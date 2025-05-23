# * `strut_twr_mountpoint::float` = [0.01,0.5,0.9], # factor of blade height where the bottom strut attaches on the tower # This puts struts at top and bottom, as a fraction of the blade position
# * `strut_bld_mountpoint::float` = [0.01,0.5,0.9], # factor of blade height where the bottom strut attaches on the blade # This puts struts at bottom 0, mid 0.5, and top 1.0 as a fraction of the blade position
# NuMad_geom_xlscsv_file_strut = nothing, or sting for the filename, or array of strings matching length of Nstruts (length(strut_twr_mountpoint))
function setupOWENS(
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

    # Unpack the tower config
    H = tower_config.H
    R = tower_config.R

    # Unpack the material config
    rhoA = material_config.rhoA
    rhoAs = material_config.rhoAs

    # Unpack the aero config
    AeroModel = aero_config.AeroModel
    AD15On = aero_config.AD15On

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

    if AddedMass_Coeff_Ca>0.0
        centrifugal_force_flag = true
    else
        centrifugal_force_flag = false
    end

    # Here is where we take the inputs from setupOWENS and break out what is going on behind the function.
    # We do some intermediate calculations on the blade shape and angles

    Nstrutperbld = length(strut_twr_mountpoint)

    if typeof(NuMad_geom_xlscsv_file_strut)==String ||
       typeof(NuMad_geom_xlscsv_file_strut)==OrderedCollections.OrderedDict{Symbol,Any}
        NuMad_geom_xlscsv_file_strut = fill(NuMad_geom_xlscsv_file_strut, Nstrutperbld)
    end

    Nbld = B
    H = maximum(shapeZ) #m,
    R = maximum(shapeX) #m,
    omega = RPM / 60 * 2 * pi
    tsr = omega*R/Vinf

    nothing

    # Here we set up the mesh using one of the pre-made meshing functions.  For this case, there is a function for the ARCUS, as 
    # well as for towered VAWTs where you can have an arbitrary blade shape with connected struts, and if the blade tips touch
    # the tower, then you can tell it to constrain them to the tower thus allowing for both H-VAWT and Darrieus designs.

    #########################################
    ### Set up mesh
    #########################################
    mesh_props = setup_mesh(mesh_config, blade_config, tower_config)


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

    # This is where the sectional properties for the tower are either read in from the file, or are directly input and could be manuplated here in the script

    #########################################
    ### Set up Sectional Properties
    #########################################

    setup_sectional_props()

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
        tower_config,
        numadIn_bld,
        numadIn_strut,
        path
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
            sections
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
            sections
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
            mass_array
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
            mass_array
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
        meshtype = meshtype,
        connectBldTips2Twr = connectBldTips2Twr,
        custommesh = custommesh,
    )
    tower_config = TowerConfig(
        H = H,
        R = R,
    )
    blade_config = BladeConfig(
        B = B,
    )
    material_config = MaterialConfig(
        rho = rho,
        mu = mu,
    )
    aero_config = AeroConfig(
        AeroModel = AeroModel,
        DynamicStallModel = DynamicStallModel,
        RPI = RPI,
        Aero_AddedMass_Active = Aero_AddedMass_Active,
        Aero_RotAccel_Active = Aero_RotAccel_Active,
        Aero_Buoyancy_Active = Aero_Buoyancy_Active,
    )

    # Run the setup function with the config structs
    
    
    return setupOWENS(
        OWENSAero,
        path;
        mesh_config = mesh_config,
        tower_config = tower_config,
        blade_config = blade_config,
        material_config = material_config,
        aero_config = aero_config,
    )
    end


function setupOWENShawt(
    OWENSAero,
    path;
    rho = 1.225,
    Nslices = 30,
    ntheta = 30,
    RPM = 1e-6,
    Vinf = 25.0,
    eta = 0.5,
    B = 3,
    H = 5.0,
    R = 2.5,
    hubR = 2.0,
    shapeZ = collect(LinRange(0, H, Nslices+1)),
    shapeX = R .* (1.0 .- 4.0 .* (shapeZ/H .- 0.5) .^ 2),#shapeX_spline(shapeZ)
    ifw = false,
    windINPfilename = "$(path)/data/turbsim/115mx115m_30x30_25.0msETM.bts",
    ifw_libfile = "$(path)/bin/libifw_c_binding",
    NuMad_geom_xlscsv_file_twr = nothing,
    NuMad_mat_xlscsv_file_twr = nothing,
    NuMad_geom_xlscsv_file_bld = nothing,
    NuMad_mat_xlscsv_file_bld = nothing,
    stack_layers_bld = nothing,
    stack_layers_scale = [1.0, 1.0],
    chord_scale = [1.0, 1.0],
    thickness_scale = [1.0, 1.0],
    Htwr_base = 2.0,
    ntelem = 10, #tower elements
    nbelem = 60, #blade elements
    ncelem = 10,
    joint_type = 2,
    c_mount_ratio = 0.05,
    AeroModel = "DMS",
    DynamicStallModel = "BV",
    RPI = true,
    biwing = false,
    hub_depth = 15.0, #Hub Beam Depth
    R_root = 10.0, # m biwing radius
    R_biwing = 30.0, # outer radius
    R_tip = 54.014, # outer radius
    nbelem_root = 30, #biwing elements, for each 
    nbelem_biwing = 30, #tip elements
    nbelem_tip = 30, #tip elements
    bshapex_root = LinRange(0.0, R_root, nbelem_root+1), #Blade shape, magnitude is relevant
    bshapez_root = zeros(nbelem_root+1), #Blade shape, magnitude is relevant
    bshapex_biwing_U = LinRange(R_root, R_biwing, nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapez_biwing_U = zeros(nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapex_biwing_L = LinRange(R_root, R_biwing, nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapez_biwing_L = zeros(nbelem_biwing+1), #Blade shape, magnitude is relevant
    bshapex_tip = LinRange(R_biwing, R_tip, nbelem_tip+1), #Blade shape, magnitude is relevant
    bshapez_tip = zeros(nbelem_tip+1), #Blade shape, magnitude is relevant
    AD15_ccw = false,
    angularOffset = 0.0,
)

    Nbld = B
    H = maximum(shapeZ) #m,
    R = maximum(shapeX) #m,
    omega = RPM / 60 * 2 * pi
    tsr = omega*R/Vinf

    # shapeX_spline = FLOWMath.Akima(shapeZ, shapeX)
    # bladelen = sum(sqrt.((shapeX[2:end].-shapeX[1:end-1]).^2 .+ (shapeZ[2:end].-shapeZ[1:end-1]).^2 ))
    # # println("bladelen 161.06: $bladelen")
    h_frac = (shapeZ[2:end] - shapeZ[1:(end-1)]) ./ shapeZ[end];
    h_elem = (shapeZ[2:end] - shapeZ[1:(end-1)])
    h = (shapeZ[2:end] + shapeZ[1:(end-1)])/2.0;

    RefArea = pi*R^2
    #########################################
    ### Set up mesh
    #########################################
    # println("Create Mesh")
    if !biwing
        mymesh, myort, myjoint, bladeIdx, bladeElem = create_hawt_mesh(;
            hub_depth = Htwr_base,
            tip_precone = H, #blade height
            R, # m bade radius
            AD15hubR = hubR,
            nblade = Nbld,
            AD15_ccw,
            ntelem, #tower elements
            nbelem, #blade elements
            bshapex = shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
            bshapez = shapeZ,
            joint_type, #hinged about y axis
            angularOffset,
        ) #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    else
        mymesh, myort, myjoint, bladeIdx, bladeElem = create_hawt_biwing_mesh(;
            nblade = Nbld,
            hub_depth,
            R_root,
            R_biwing,
            R_tip,
            ntelem,
            nbelem_root,
            nbelem_biwing,
            nbelem_tip,
            bshapex_root,
            bshapez_root,
            bshapex_biwing_U,
            bshapez_biwing_U,
            bshapex_biwing_L,
            bshapez_biwing_L,
            bshapex_tip,
            bshapez_tip,
            joint_type,
            angularOffset = -pi/2,
        )
    end

    nTwrElem = Int(mymesh.meshSeg[1])+1

    # PyPlot.figure()
    # PyPlot.plot(mymesh.x,mymesh.z,"b-")
    #  for myi = 1:length(mymesh.x)
    #      PyPlot.text(mymesh.x[myi].+rand()/30,mymesh.z[myi].+rand()/30,"$myi",ha="center",va="center")
    #      PyPlot.draw()
    #      sleep(0.1)
    #  end
    # PyPlot.xlabel("x")
    # PyPlot.ylabel("y")
    # PyPlot.axis("equal")

    #########################################
    ### Set up Sectional Properties
    #########################################
    # println("Calculate/Set up sectional properties")
    #Tower
    if !isnothing(NuMad_geom_xlscsv_file_twr)
        numadIn_twr = readNuMadGeomCSV(NuMad_geom_xlscsv_file_twr)
    else
        n_web = 0
        n_stack = 2
        n_segments = 2
        span = [
            0.0,
            6.607421057,
            13.21484211,
            19.82226317,
            26.42968423,
            33.03710529,
            39.64452634,
            46.2519474,
            52.85936846,
            59.46678951,
            66.07421057,
            72.68163163,
            79.28905268,
            85.89647374,
            92.5038948,
            99.11131586,
            105.7187369,
            112.326158,
            118.933579,
            125.5410001,
            132.1484211,
        ]
        airfoil = [
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
            "circular",
        ]
        te_type = [
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
            "round",
        ]
        twist_d = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
        chord = [
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            0.25,
            0.25,
            0.25,
        ]
        xoffset = [
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
        ]
        aerocenter = [
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
            0.25,
        ]
        stack_mat_types = [8, 2]
        stack_layers = [
            70 3;
            70 3;
            70 3;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03;
            30.303 3.03
        ]
        segments = [
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0;
            -1.0 0.0 1.0
        ]
        DPtypes = [
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" "";
            "" "" ""
        ]
        skin_seq = [
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2]);
            Seq([2, 1, 2]) Seq([2, 1, 2])
        ]
        web_seq = Array{Seq,2}(undef, length(twist_d), 0) #can be any number of stack nums, so we have to make non-square containers
        web_dp = Array{Seq,2}(undef, length(twist_d), 0) #this is fixed size square, but it's easier to do it this way

        numadIn_twr = NuMad(
            n_web,
            n_stack,
            n_segments,
            span,
            airfoil,
            te_type,
            twist_d,
            chord,
            xoffset,
            aerocenter,
            stack_mat_types,
            stack_layers,
            segments,
            DPtypes,
            skin_seq,
            web_seq,
            web_dp,
        )
    end

    #Add the full path
    for (i, airfoil) in enumerate(numadIn_twr.airfoil)
        numadIn_twr.airfoil[i] = "$path/Airfoils/$airfoil"
    end

    if !isnothing(NuMad_mat_xlscsv_file_twr)
        plyprops_twr = readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_twr)
    else
        names = [
            "CLA_5500",
            "CBX_2400",
            "ETLX_2400",
            "Airex_C70_55",
            "EBX_2400_x10",
            "ETLX_2400_x10",
            "Airex_C70_55_x10",
            "CFP-baseline",
        ]
        plies = [
            Composites.Material{Float64}(
                9.824e10,
                5.102e9,
                4.274e9,
                0.3,
                1540.0,
                8.75634139e8,
                5.92949102e8,
                1.0e8,
                1.0e8,
                1.0e8,
                0.00066,
            ),
            Composites.Material{Float64}(
                1.4931e10,
                1.4931e10,
                2.389e10,
                0.3,
                1530.0,
                4.55053962e8,
                4.55053962e8,
                1.0e8,
                1.0e8,
                1.0e8,
                0.0008100000000000001,
            ),
            Composites.Material{Float64}(
                2.0333e10,
                9.305e9,
                4.756e9,
                0.3,
                1900.0,
                5.30896289e8,
                5.30896289e8,
                1.0e8,
                1.0e8,
                1.0e8,
                0.00066,
            ),
            Composites.Material{Float64}(
                4.5e7,
                4.5e7,
                2.2e7,
                0.2,
                59.0,
                1.0e8,
                1.0e8,
                1.0e8,
                1.0e8,
                1.0e8,
                0.001,
            ),
            Composites.Material{Float64}(
                9.824e11,
                5.102e10,
                4.274e10,
                0.3,
                15300.0,
                4.55053962e9,
                4.55053962e9,
                1.0e8,
                1.0e8,
                1.0e8,
                7.000000000000001e-5,
            ),
            Composites.Material{Float64}(
                1.4931e11,
                1.4931e11,
                2.389e11,
                0.3,
                19000.0,
                5.30896289e9,
                5.30896289e9,
                1.0e8,
                1.0e8,
                1.0e8,
                8.0e-5,
            ),
            Composites.Material{Float64}(
                2.03335e11,
                9.3051e10,
                4.756e10,
                0.2,
                590.0,
                1.0e9,
                1.0e9,
                1.0e8,
                1.0e8,
                1.0e8,
                7.000000000000001e-5,
            ),
            Composites.Material{Float64}(
                1.576e11,
                9.1e9,
                3.3e9,
                0.263,
                1600.0,
                2.236e9,
                1.528e9,
                1.0e8,
                1.0e8,
                1.0e8,
                0.00066,
            ),
        ]
        plyprops_twr = plyproperties(names, plies)
    end

    twr_precompoutput, twr_precompinput, lam_U_twr, lam_L_twr, lam_W_twr =
        getOWENSPreCompOutput(numadIn_twr; plyprops = plyprops_twr)
    sectionPropsArray_twr = getSectPropsFromOWENSPreComp(
        LinRange(0, 1, nTwrElem),
        numadIn_twr,
        twr_precompoutput;
        precompinputs = twr_precompinput,
    )
    stiff_twr, mass_twr = getSectPropsFromOWENSPreComp(
        LinRange(0, 1, nTwrElem),
        numadIn_twr,
        twr_precompoutput;
        GX = true,
    )

    #Blades
    if !isnothing(NuMad_geom_xlscsv_file_bld)
        numadIn_bld = readNuMadGeomCSV(NuMad_geom_xlscsv_file_bld)
    else
        n_web = 1
        n_stack = 7
        n_segments = 12
        span = [
            0.0,
            6.607,
            13.215,
            19.822,
            26.43,
            33.037,
            39.645,
            46.252,
            52.859,
            59.467,
            66.074,
            72.682,
            79.289,
            85.896,
            92.504,
            99.111,
            105.719,
            112.326,
            118.934,
            125.541,
            132.148,
        ]
        airfoil = [
            "Cylinder1",
            "Cylinder1",
            "Cylinder1",
            "Cylinder2",
            "Cylinder2",
            "Cylinder2",
            "DU21_A17",
            "DU21_A17",
            "DU21_A17",
            "DU25_A17",
            "DU25_A17",
            "DU25_A17",
            "DU30_A17",
            "DU30_A17",
            "DU30_A17",
            "DU35_A17",
            "DU35_A17",
            "DU35_A17",
            "DU40_A17",
            "DU40_A17",
            "NACA64_A17",
        ]
        te_type = [
            "round",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
            "sharp",
        ]
        twist_d = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
        chord = [
            10.0,
            10.0,
            9.0,
            8.0,
            8.0,
            7.0,
            7.0,
            6.0,
            6.0,
            6.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
        ]
        xoffset = [
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
            0.3,
        ]
        aerocenter = [
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
        ]
        stack_mat_types = [8, 2, 4, 8, 8, 8, 4]
        if isnothing(stack_layers_bld)
            stack_layers = [
                30.0 2.0 15.0 25.0 25.0 2.0 13.0;
                15.0 2.0 10.0 13.0 11.0 2.0 11.0;
                10.0 1.0 8.0 10.0 10.0 2.0 10.0;
                8.0 1.0 6.0 9.0 10.0 1.0 9.0;
                7.0 1.0 5.0 8.0 9.0 1.0 7.0;
                6.0 1.0 4.0 8.0 9.0 1.0 6.0;
                6.0 1.0 4.0 8.0 8.0 1.0 5.0;
                6.0 1.0 4.0 7.0 7.0 1.0 5.0;
                7.0 1.0 3.0 6.0 6.0 1.0 5.0;
                8.0 1.0 3.0 6.0 6.0 1.0 5.0;
                8.0 1.0 3.0 6.0 6.0 1.0 5.0;
                7.0 1.0 3.0 6.0 6.0 1.0 5.0;
                7.0 1.0 3.0 6.0 6.0 2.0 5.0;
                7.0 1.0 3.0 6.0 6.0 2.0 5.0;
                7.0 1.0 3.0 7.0 8.0 3.0 5.0;
                7.0 2.0 3.0 9.0 12.0 3.0 6.0;
                10.0 3.0 4.0 11.0 15.0 3.0 6.0;
                12.0 3.0 4.0 13.0 15.0 3.0 6.0;
                12.0 3.0 4.0 15.0 15.0 3.0 6.0;
                12.0 3.0 4.0 15.0 15.0 3.0 6.0;
                10.0 1.0 4.0 10.0 12.0 1.0 5.0
            ]
        else
            stack_layers = stack_layers_bld
        end
        segments = [
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0;
            -1.0 -0.95 -0.5 -0.3 -0.1 -0.095 0.0 0.095 0.1 0.3 0.5 0.95 1.0
        ]
        DPtypes = [
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" "";
            "" "" "" "" "" "" "" "" "" "" "" "" ""
        ]
        skin_seq = [
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2]);
            Seq([2, 5, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 1, 2]) Seq([2, 1, 2]) Seq([2, 3, 2]) Seq([2, 4, 2]) Seq([2, 4, 2]) Seq([2, 3, 2]) Seq([2, 5, 2])
        ]
        web_seq = [
            Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]); Seq([6, 7, 6]);;
        ]
        web_dp = [
            Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]); Seq([9, 3, 3, 9]);;
        ]

        numadIn_bld = NuMad(
            n_web,
            n_stack,
            n_segments,
            span,
            airfoil,
            te_type,
            twist_d,
            chord,
            xoffset,
            aerocenter,
            stack_mat_types,
            stack_layers,
            segments,
            DPtypes,
            skin_seq,
            web_seq,
            web_dp,
        )
    end
    for icol = 1:length(numadIn_bld.stack_layers[1, :])
        numadIn_bld.stack_layers[:, icol] .*= LinRange(
            stack_layers_scale[1],
            stack_layers_scale[2],
            length(numadIn_bld.chord),
        )
    end
    numadIn_bld.chord .*=
        LinRange(chord_scale[1], chord_scale[2], length(numadIn_bld.chord))

    for (i, airfoil) in enumerate(numadIn_bld.airfoil)
        numadIn_bld.airfoil[i] = "$path/Airfoils/$airfoil"
    end

    if !isnothing(NuMad_mat_xlscsv_file_bld)
        plyprops_bld = readNuMadMaterialsCSV(NuMad_mat_xlscsv_file_bld)
    else
        names = [
            "CLA_5500",
            "CBX_2400",
            "ETLX_2400",
            "Airex_C70_55",
            "EBX_2400_x10",
            "ETLX_2400_x10",
            "Airex_C70_55_x10",
            "CFP-baseline",
        ]
        plies = [
            Composites.Material{Float64}(
                9.824e10,
                5.102e9,
                4.274e9,
                0.3,
                1540.0,
                8.75634139e8,
                5.92949102e8,
                1.0e8,
                1.0e8,
                1.0e8,
                0.00066,
            ),
            Composites.Material{Float64}(
                1.4931e10,
                1.4931e10,
                2.389e10,
                0.3,
                1530.0,
                4.55053962e8,
                4.55053962e8,
                1.0e8,
                1.0e8,
                1.0e8,
                0.0008100000000000001,
            ),
            Composites.Material{Float64}(
                2.0333e10,
                9.305e9,
                4.756e9,
                0.3,
                1900.0,
                5.30896289e8,
                5.30896289e8,
                1.0e8,
                1.0e8,
                1.0e8,
                0.00066,
            ),
            Composites.Material{Float64}(
                4.5e7,
                4.5e7,
                2.2e7,
                0.2,
                59.0,
                1.0e8,
                1.0e8,
                1.0e8,
                1.0e8,
                1.0e8,
                0.001,
            ),
            Composites.Material{Float64}(
                9.824e11,
                5.102e10,
                4.274e10,
                0.3,
                15300.0,
                4.55053962e9,
                4.55053962e9,
                1.0e8,
                1.0e8,
                1.0e8,
                7.000000000000001e-5,
            ),
            Composites.Material{Float64}(
                1.4931e11,
                1.4931e11,
                2.389e11,
                0.3,
                19000.0,
                5.30896289e9,
                5.30896289e9,
                1.0e8,
                1.0e8,
                1.0e8,
                8.0e-5,
            ),
            Composites.Material{Float64}(
                2.03335e11,
                9.3051e10,
                4.756e10,
                0.2,
                590.0,
                1.0e8,
                1.0e8,
                1.0e8,
                1.0e8,
                1.0e8,
                7.000000000000001e-5,
            ),
            Composites.Material{Float64}(
                1.576e11,
                9.1e9,
                3.3e9,
                0.263,
                1600.0,
                2.236e9,
                1.528e9,
                1.0e8,
                1.0e8,
                1.0e8,
                0.00066,
            ),
        ]
        plyprops_bld = plyproperties(names, plies)
    end

    # Get blade spanwise position
    bld1start = Int(mymesh.structuralNodeNumbers[1, 1])
    bld1end = Int(mymesh.structuralNodeNumbers[1, end])
    spanpos = [
        0.0;
        cumsum(
            sqrt.(
                diff(mymesh.y[bld1start:bld1end]) .^ 2 .+
                diff(mymesh.z[bld1start:bld1end]) .^ 2,
            ),
        )
    ]

    #TODO: fix this, it is a hack
    if biwing
        spanpos = collect(LinRange(0, spanpos[end], length(spanpos)-3))
    end

    if length(thickness_scale)==2
        yscale = collect(
            LinRange(thickness_scale[1], thickness_scale[2], length(numadIn_bld.span)),
        )
    elseif length(thickness_scale)==length(numadIn_bld.span)
        yscale = thickness_scale
    end

    bld_precompoutput, bld_precompinput, lam_U_bld, lam_L_bld, lam_W_bld =
        getOWENSPreCompOutput(numadIn_bld; yscale, plyprops = plyprops_bld)
    sectionPropsArray_bld = getSectPropsFromOWENSPreComp(
        spanpos,
        numadIn_bld,
        bld_precompoutput;
        precompinputs = bld_precompinput,
    )
    stiff_bld, mass_bld =
        getSectPropsFromOWENSPreComp(spanpos, numadIn_bld, bld_precompoutput; GX = true)

    #Struts
    # They are the same as the end properties of the blades

    # Combined Section Props
    bldssecprops = collect(Iterators.flatten(fill(sectionPropsArray_bld, Nbld)))
    sectionPropsArray =
        [fill(sectionPropsArray_twr[1], length(sectionPropsArray_twr)); bldssecprops]#;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str;sectionPropsArray_str]
    # sectionPropsArray = fill(sectionPropsArray[end],length(sectionPropsArray))
    rotationalEffects = ones(mymesh.numEl)

    # GXBeam sectional properties
    stiff_blds = collect(Iterators.flatten(fill(stiff_bld, Nbld)))
    stiff_array = [stiff_twr; stiff_blds]

    mass_blds = collect(Iterators.flatten(fill(mass_bld, Nbld)))
    mass_array = [mass_twr; mass_blds]

    system, assembly, sections = owens_to_gx(
        mymesh,
        myort,
        myjoint,
        sectionPropsArray,
        stiff_array,
        mass_array;
        VTKmeshfilename = "$path/vtk/HAWT_biwing",
    )

    #store data in element object
    myel = OWENSFEA.El(
        sectionPropsArray,
        myort.Length,
        myort.Psi_d,
        myort.Theta_d,
        myort.Twist_d,
        rotationalEffects,
    )

    #########################################
    ### Set up aero forces
    #########################################
    # println("Initialize Aerodynamics")
    # chord_spl = safeakima(numadIn_bld.span./maximum(numadIn_bld.span), numadIn_bld.chord,LinRange(0,1,Nslices))
    # OWENSAero.setupTurb(shapeX,shapeZ,B,chord_spl,tsr,Vinf;AeroModel,DynamicStallModel,
    # afname = "$path/Airfoils/NACA_0021.dat", #TODO: map to the numad input
    # ifw,
    # windINPfilename,
    # ifw_libfile,
    # ntheta,Nslices,rho,eta,RPI)

    aeroForces(t, azi) = mapACDMS(
        t,
        azi,
        mymesh,
        myel,
        OWENSAero.AdvanceTurbineInterpolate;
        alwaysrecalc = true,
    )


    # Calculate mass breakout of each material
    mass_breakout_bld = get_material_mass(plyprops_bld, numadIn_bld)
    mass_breakout_blds = mass_breakout_bld .* length(mymesh.structuralNodeNumbers[:, 1])
    mass_breakout_twr =
        get_material_mass(plyprops_twr, numadIn_twr; int_start = 0.0, int_stop = Htwr_base)

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
    aeroForces,
    mass_breakout_blds,
    mass_breakout_twr,
    bladeIdx,
    bladeElem,
    system,
    assembly,
    sections
end
