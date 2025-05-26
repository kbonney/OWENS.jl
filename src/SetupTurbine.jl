using ..OWENS
using ..OWENSFEA
using ..OWENSOpenFASTWrappers
using ..OWENSAero
using ..Composites
using ..FLOWMath
using ..QuadGK
using ..OrderedCollections

# Import types and functions from setup_utilities.jl
import ..OWENS: MeshConfig, TowerConfig, BladeConfig, MaterialConfig, AeroConfig
import ..OWENS: default_mesh_config, default_tower_config, default_blade_config, default_material_config, default_aero_config
import ..OWENS: setup_mesh, setup_sectional_props, setup_aerodynamic_model

# * `strut_twr_mountpoint::float` = [0.01,0.5,0.9], # factor of blade height where the bottom strut attaches on the tower # This puts struts at top and bottom, as a fraction of the blade position
# * `strut_bld_mountpoint::float` = [0.01,0.5,0.9], # factor of blade height where the bottom strut attaches on the blade # This puts struts at bottom 0, mid 0.5, and top 1.0 as a fraction of the blade position
# NuMad_geom_xlscsv_file_strut = nothing, or sting for the filename, or array of strings matching length of Nstruts (length(strut_twr_mountpoint))
function setupOWENS_struct(
    OWENSAero::Module,
    path::String;
    verbosity = 1,
    VTKmeshfilename = nothing,
    return_componentized = false,
    mesh_config = default_mesh_config(),
    tower_config = default_tower_config(),
    blade_config = default_blade_config(),
    material_config = default_material_config(),
    aero_config = default_aero_config(),
)

    # Unpack the mesh config
    meshtype = mesh_config.meshtype
    connectBldTips2Twr = mesh_config.connectBldTips2Twr

    # Unpack the blade config
    B = blade_config.B
    shapeX = blade_config.shapeX
    shapeZ = blade_config.shapeZ
    shapeY = blade_config.shapeY

    # Unpack the tower config

    # Unpack the material config


    # Unpack the aero config
    AeroModel = aero_config.AeroModel

    custom_mesh_outputs = []

    if AeroModel=="AD"
        AD15On = true
    else
        AD15On = false
    end

    if meshtype == "Darrieus"
        connectBldTips2Twr = true
    else
        connectBldTips2Twr = false
    end

    if minimum(shapeZ)!=0
        @error "blade shapeZ must start at 0.0"
    end

    if material_config.AddedMass_Coeff_Ca>0.0
        centrifugal_force_flag = true
    else
        centrifugal_force_flag = false
    end

    # Here is where we take the inputs from setupOWENS and break out what is going on behind the function.
    # We do some intermediate calculations on the blade shape and angles

    Nstrutperbld = length(tower_config.strut_twr_mountpoint)

    NuMad_geom_xlscsv_file_strut = tower_config.NuMad_geom_xlscsv_file_strut
    NuMad_mat_xlscsv_file_strut = tower_config.NuMad_mat_xlscsv_file_strut

    if typeof(NuMad_geom_xlscsv_file_strut)==String ||
       typeof(NuMad_geom_xlscsv_file_strut)==OrderedCollections.OrderedDict{Symbol,Any}
        NuMad_geom_xlscsv_file_strut = fill(NuMad_geom_xlscsv_file_strut, Nstrutperbld)
    end

    Nbld = B
    H = maximum(shapeZ) #m,
    R = 0.0 # Set radius to 0.0 for default case
    omega = aero_config.RPM / 60 * 2 * pi
    tsr = omega*R/aero_config.Vinf

    nothing

    # Here we set up the mesh using one of the pre-made meshing functions.  For this case, there is a function for the ARCUS, as 
    # well as for towered VAWTs where you can have an arbitrary blade shape with connected struts, and if the blade tips touch
    # the tower, then you can tell it to constrain them to the tower thus allowing for both H-VAWT and Darrieus designs.

    #########################################
    ### Set up mesh
    #########################################
    mesh_props = setup_mesh(mesh_config, blade_config, tower_config, verbosity)

    mymesh = mesh_props.mymesh
    myort = mesh_props.myort
    myjoint = mesh_props.myjoint

    nTwrElem = Int(mymesh.meshSeg[1])+1
    # try
    #     if contains(NuMad_mat_xlscsv_file_bld,"34m") || meshtype == "ARCUS" #TODO: this is really odd, 
    #         nTwrElem = Int(mymesh.meshSeg[1])+1
    #     end
    # catch
    #     nTwrElem = Int(mymesh.meshSeg[1])
    # end

    nothing

    # Here is a way that you can visualize the nodal numbers of the mesh

    # PyPlot.figure()
    # for icon = 1:length(mymesh.conn[:,1])
    #     idx1 = mymesh.conn[icon,1]
    #     idx2 = mymesh.conn[icon,2]
    #     PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"k.-")
    #     PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
    #     # sleep(0.1)
    # end

    # for ijoint = 1:length(myjoint[:,1])
    #     idx2 = Int(myjoint[ijoint,2])
    #     idx1 = Int(myjoint[ijoint,3])
    #     PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"r.-")
    #     # sleep(0.1)
    # end

    # Unpack to get the element numbers for the extra components
    components,
    s2t_idx,
    intra_blade_cable_startidx_el,
    intra_blade_cable_endidx_el,
    topel_idx,
    s2b_idx = mesh_props.custom_mesh_outputs

    # Unpack AD15bldNdIdxRng from mesh_props
    AD15bldNdIdxRng = mesh_props.AD15bldNdIdxRng
    AD15bldElIdxRng = mesh_props.AD15bldElIdxRng

    # This is where the sectional properties for the tower are either read in from the file, or are directly input and could be manuplated here in the script

    #########################################
    ### Set up Sectional Properties
    #########################################

    sectional_props = setup_sectional_props(
        components,
        mymesh,
        path,
        aero_config.rho,
        material_config.AddedMass_Coeff_Ca,
        tower_config.NuMad_geom_xlscsv_file_twr,
        tower_config.NuMad_mat_xlscsv_file_twr,
        blade_config.NuMad_geom_xlscsv_file_bld,
        blade_config.NuMad_mat_xlscsv_file_bld,
        tower_config.NuMad_geom_xlscsv_file_strut,
        tower_config.NuMad_mat_xlscsv_file_strut,
        nothing, # NuMad_geom_xlscsv_file_intra_blade_cable
        nothing, # NuMad_mat_xlscsv_file_intra_blade_cable
        nothing, # NuMad_geom_xlscsv_file_guys
        nothing, # NuMad_mat_xlscsv_file_guys
        verbosity
    )
    sectionPropsArray = sectional_props.sectionPropsArray
    stiff_array = sectional_props.stiff_array
    mass_array = sectional_props.mass_array
    rotationalEffects = sectional_props.rotationalEffects

    if verbosity>1
        println(
            "\nTotal Mass: $(mapreduce(sum, +, [components[icomponent].mass for icomponent = 1:size(components)[1]])) kg",
        )
        println(
            "Total Material Cost: \$$(mapreduce(sum, +, [components[icomponent].mass for icomponent = 1:size(components)[1]]))",
        )
    end

    # # Plot the rotational effects' status
    # for iel = 1:length(mymesh.conn[:,1])
    #     if rotationalEffects[iel] == 1
    #         idx1 = mymesh.conn[iel,1]
    #         idx2 = mymesh.conn[iel,2]
    #         PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],linewidth=10.0,"b.-")
    #         # PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
    #         # sleep(0.1)
    #     end
    # end

    if length(sectionPropsArray)<mymesh.numEl
        @warn "There are more mesh elements $(mymesh.numEl) than sectional properties $(length(sectionPropsArray)), applying the last element's sectional properties to the remaining"
        n_diff = mymesh.numEl - length(sectionPropsArray)
        sectionPropsArray = [sectionPropsArray; fill(sectionPropsArray[end], n_diff)]
        stiff_array = [stiff_array; fill(stiff_array[end], n_diff)]
        mass_array = [mass_array; fill(mass_array[end], n_diff)]
    end

    #### store data in element object
    myel = OWENSFEA.El(
        sectionPropsArray,
        myort.Length,
        myort.Psi_d,
        myort.Theta_d,
        myort.Twist_d,
        rotationalEffects,
    )

    system, assembly, sections = OWENS.owens_to_gx(
        mymesh,
        myort,
        myjoint,
        sectionPropsArray,
        stiff_array,
        mass_array;
        VTKmeshfilename,
    )

    function calcMass(sectionPropsArray, myort)
        mass = 0.0
        lens = zeros(length(sectionPropsArray))
        rhoAs = zeros(length(sectionPropsArray))
        for (i, sectionProp) in enumerate(sectionPropsArray)
            lenEl = 0
            try
                lenEl = myort.Length[i]
            catch
                lenEl = myort.Length[i-1]
            end
            lens[i] = lenEl
            rhoA = sectionProp.rhoA[1]
            rhoAs[i] = rhoA
            mass += lenEl*rhoA
        end
        span_array = cumsum(lens)
        spl_rhoa = FLOWMath.Akima(span_array, rhoAs)
        mass2, error = QuadGK.quadgk(spl_rhoa, span_array[1], span_array[end], atol = 1e-10)
        return mass, mass2
    end

    turb_masskg1, turb_masskg = calcMass(myel.props, myort)
    if verbosity>1
        println("\nMass of Turbine: $turb_masskg kg, $turb_masskg1 kg")
    end

    # Set up the OWENSAero aerodynamics if used
    blade_component_index = findfirst(s -> s.name[1:(end-1)] == "blade", components) #This assumes 9 or fewer blades, not including struts
    numadIn_bld = components[blade_component_index].nuMadIn

    strut_component_index = findall(s -> contains(s.name, "strut"), components) #This assumes 9 or fewer blades, not including struts
    numadIn_strut = [
        components[strut_component_index1].nuMadIn for
        strut_component_index1 in strut_component_index[1:Nbld:end]
    ]
    
    aeroForcesAD, deformAeroAD, aeroForcesACDMS, deformAeroACDMS = setup_aerodynamic_model(
        blade_config,
        aero_config,
        # tower_config,
        mesh_config,
        mesh_props,
        components,
        numadIn_bld,
        numadIn_strut,
        path,
        verbosity,
    )

    if return_componentized
        if AD15On
            return mymesh,
            myel,
            myort,
            myjoint,
            components,
            aeroForcesAD,
            deformAeroAD,
            system,
            assembly,
            sections,
            R
        else
            return mymesh,
            myel,
            myort,
            myjoint,
            components,
            aeroForcesACDMS,
            deformAeroACDMS,
            system,
            assembly,
            sections,
            R
        end
    else
        @warn "Not using the componetized model is being depreciated, please consider updating to return_componentized=true"

        mass_twr = []
        stiff_twr = []
        twr_precompinput = []
        twr_precompoutput = []
        plyprops_twr = []
        numadIn_twr = []
        lam_U_twr = []
        lam_L_twr = []
        mass_breakout_twr = []
        mass_bld = []
        stiff_bld = []
        bld_precompinput = []
        bld_precompoutput = []
        plyprops_bld = []
        numadIn_bld = []
        lam_U_bld = []
        lam_L_bld = []
        mass_breakout_blds = []

        for component in components
            if contains(component.name, "tower")
                mass_twr = component.mass_matrix
                stiff_twr = component.stiff_matrix
                twr_precompinput = component.preCompInput
                twr_precompoutput = component.preCompOutput
                plyprops_twr = component.plyProps
                numadIn_twr = component.nuMadIn
                lam_U_twr = component.lam_U
                lam_L_twr = component.lam_L
                mass_breakout_twr = component.mass
            end

            if contains(component.name, "blade1")
                mass_bld = component.mass_matrix
                stiff_bld = component.stiff_matrix
                bld_precompinput = component.preCompInput
                bld_precompoutput = component.preCompOutput
                plyprops_bld = component.plyProps
                numadIn_bld = component.nuMadIn
                lam_U_bld = component.lam_U
                lam_L_bld = component.lam_L
                mass_breakout_blds = component.mass
            end
        end

        if AD15On
            return mymesh,
            myel,
            myort,
            myjoint,
            sectionPropsArray,
            mass_twr,
            mass_bld,
            stiff_twr,
            stiff_bld,
            bld_precompinput,
            bld_precompoutput,
            plyprops_bld,
            numadIn_bld,
            lam_U_bld,
            lam_L_bld,
            twr_precompinput,
            twr_precompoutput,
            plyprops_twr,
            numadIn_twr,
            lam_U_twr,
            lam_L_twr,
            aeroForcesAD,
            deformAeroAD,
            mass_breakout_blds,
            mass_breakout_twr,
            system,
            assembly,
            sections,
            AD15bldNdIdxRng,
            AD15bldElIdxRng,
            custom_mesh_outputs,
            stiff_array,
            mass_array,
            R
        else
            return mymesh,
            myel,
            myort,
            myjoint,
            sectionPropsArray,
            mass_twr,
            mass_bld,
            stiff_twr,
            stiff_bld,
            bld_precompinput,
            bld_precompoutput,
            plyprops_bld,
            numadIn_bld,
            lam_U_bld,
            lam_L_bld,
            twr_precompinput,
            twr_precompoutput,
            plyprops_twr,
            numadIn_twr,
            lam_U_twr,
            lam_L_twr,
            aeroForcesACDMS,
            deformAeroACDMS,
            mass_breakout_blds,
            mass_breakout_twr,
            system,
            assembly,
            sections,
            AD15bldNdIdxRng,
            AD15bldElIdxRng,
            custom_mesh_outputs,
            stiff_array,
            mass_array,
            R
        end
    end
end

function setupOWENS(
    OWENSAero,
    path;
    rho = 1.225,
    mu = 1.7894e-5,
    Nslices = 30,
    ntheta = 30,
    RPM = 1e-6,
    Vinf = 25.0,
    eta = 0.5,
    B = 3,
    H = 5.0,
    R = 2.5,
    shapeZ = collect(LinRange(0, H, Nslices+1)),
    shapeX = R .* (1.0 .- 4.0 .* (shapeZ/H .- 0.5) .^ 2),#shapeX_spline(shapeZ)
    ntelem = 10, #tower elements
    nbelem = 60, #blade elements
    ncelem = 10,
    nselem = 5,
    shapeY = zeros(Nslices+1),
    ifw = false,
    AD15hubR = 0.1,
    WindType = 1,
    delta_t = 0.01,
    numTS = 100,
    adi_lib = "$(path)../../../../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding",
    adi_rootname = "./Example",
    windINPfilename = "$(path)/data/turbsim/115mx115m_30x30_25.0msETM.bts",
    ifw_libfile = "$(path)/bin/libifw_c_binding",
    NuMad_geom_xlscsv_file_twr = nothing,
    NuMad_mat_xlscsv_file_twr = nothing,
    NuMad_geom_xlscsv_file_bld = nothing,
    NuMad_mat_xlscsv_file_bld = nothing,
    NuMad_geom_xlscsv_file_strut = nothing,
    NuMad_mat_xlscsv_file_strut = nothing,
    NuMad_geom_xlscsv_file_intra_blade_cable = nothing,
    NuMad_mat_xlscsv_file_intra_blade_cable = nothing,
    NuMad_geom_xlscsv_file_guys = nothing,
    NuMad_mat_xlscsv_file_guys = nothing,
    stack_layers_bld = nothing,
    stack_layers_scale = [1.0, 1.0],
    chord_scale = [1.0, 1.0],
    thickness_scale = [1.0, 1.0],
    Htwr_base = 2.0,
    Htwr_blds = H,
    strut_twr_mountpoint = [0.25, 0.75],
    strut_bld_mountpoint = [0.25, 0.75],
    joint_type = 2,
    c_mount_ratio = 0.05,
    angularOffset = -pi/2,
    AeroModel = "DMS",
    DynamicStallModel = "BV",
    RPI = true,
    Aero_AddedMass_Active = false,
    Aero_RotAccel_Active = false,
    Aero_Buoyancy_Active = false,
    cables_connected_to_blade_base = true,
    meshtype = "Darrieus",
    custommesh = nothing,
    AddedMass_Coeff_Ca = 0.0,
    verbosity = 1,
    VTKmeshfilename = nothing,
    return_componentized = false,
)

    # Bundle the inputs into config structs
    mesh_config = MeshConfig(
        Nslices = Nslices,
        ntheta = ntheta,
        ntelem = ntelem,
        nbelem = nbelem,
        ncelem = ncelem,
        nselem = nselem,
        meshtype = meshtype,
        custommesh = custommesh,
        connectBldTips2Twr = meshtype == "Darrieus", # Connect blade tips to tower for Darrieus type
        AD15_ccw = true,
    )

    tower_config = TowerConfig(
        Htwr_base = Htwr_base,
        Htwr_blds = Htwr_blds,
        strut_twr_mountpoint = strut_twr_mountpoint,
        strut_bld_mountpoint = strut_bld_mountpoint,
        joint_type = joint_type,
        c_mount_ratio = c_mount_ratio,
        angularOffset = angularOffset,
        NuMad_geom_xlscsv_file_twr = NuMad_geom_xlscsv_file_twr,
        NuMad_mat_xlscsv_file_twr = NuMad_mat_xlscsv_file_twr,
    )

    blade_config = BladeConfig(
        B = B,
        H = H,
        R = R,
        shapeZ = shapeZ,
        shapeX = shapeX,
        shapeY = shapeY,
        NuMad_geom_xlscsv_file_bld = NuMad_geom_xlscsv_file_bld,
        NuMad_mat_xlscsv_file_bld = NuMad_mat_xlscsv_file_bld,
    )

    material_config = MaterialConfig(
        stack_layers_bld = stack_layers_bld,
        stack_layers_scale = stack_layers_scale,
        chord_scale = chord_scale,
        thickness_scale = thickness_scale,
        AddedMass_Coeff_Ca = AddedMass_Coeff_Ca,
    )

    aero_config = AeroConfig(
        rho = rho,
        mu = mu,
        RPM = RPM,
        Vinf = Vinf,
        eta = eta,
        delta_t = delta_t,
        AD15hubR = AD15hubR,
        WindType = WindType,
        AeroModel = AeroModel,
        DynamicStallModel = DynamicStallModel,
        numTS = numTS,
        adi_lib = adi_lib,
        adi_rootname = adi_rootname,
        windINPfilename = windINPfilename,
        ifw_libfile = ifw_libfile,
        ifw = ifw,
        RPI = RPI,
        Aero_AddedMass_Active = Aero_AddedMass_Active,
        Aero_RotAccel_Active = Aero_RotAccel_Active,
        Aero_Buoyancy_Active = Aero_Buoyancy_Active,
        centrifugal_force_flag = AddedMass_Coeff_Ca > 0.0,
    )

    # Run the setup function with the config structs
    return setupOWENS_struct(
        OWENSAero,
        path;
        verbosity = verbosity,
        VTKmeshfilename = VTKmeshfilename,
        return_componentized = return_componentized,
        mesh_config = mesh_config,
        tower_config = tower_config,
        blade_config = blade_config,
        material_config = material_config,
        aero_config = aero_config,
    )
end
