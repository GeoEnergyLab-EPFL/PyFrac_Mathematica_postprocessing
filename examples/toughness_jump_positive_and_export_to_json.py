# -*- coding: utf-8 -*-
"""
This file is part of PyFrac.

Created by Carlo Peruzzo on Fri Apr 17 23:16:25 2020.
Copyright (c) "ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory", 2016-2020.
All rights reserved. See the LICENSE.TXT file for more details.
"""

# local imports
from mesh import CartesianMesh
from properties import MaterialProperties, FluidProperties, InjectionProperties, SimulationProperties
from fracture import Fracture
from controller import Controller
from fracture_initialization import Geometry, InitializationParameters
from utility import setup_logging_to_console

# setting up the verbosity level of the log at console
setup_logging_to_console(verbosity_level='debug')

# creating mesh
Mesh = CartesianMesh(0.25, 0.15, 69, 39)

# solid properties
nu = 0.4                            # Poisson's ratio
youngs_mod = 3.3e10                 # Young's modulus
Eprime = youngs_mod / (1 - nu ** 2) # plain strain modulus
K_Ic1 = 5.6e6                       # fracture toughness


def My_KIc_func(x, y):
    """ The function providing the fracture toughness"""
    if x<0.2:
        return K_Ic1
    else:
        return K_Ic1*3.6

Solid = MaterialProperties(Mesh,
                           Eprime,
                           K1c_func=My_KIc_func,
                           minimum_width=1e-8)

# injection parameters
Q0 = 0.05  # injection rate
Injection = InjectionProperties(Q0, Mesh, source_coordinates=[0.1,0.])

# fluid properties
Fluid = FluidProperties(viscosity=1.1e-3)

# simulation properties
simulProp = SimulationProperties()
simulProp.finalTime = 0.0003                                               # the time at which the simulation stops
myfolder="./Data/toughness_jump_3p6"
simulProp.set_outputFolder(myfolder) # the disk address where the files are saved
simulProp.projMethod = 'LS_continousfront'
simulProp.frontAdvancing = 'implicit'
simulProp.useBlockToeplizCompression = True
simulProp.saveToDisk = False
simulProp.bckColor = 'K1c'
simulProp.saveFluidVelAsVector = True

# To decide what you will see when you print:
#simulProp.plotVar = ['ffvf','regime']
#simulProp.plotVar = ['footprint','regime']
#simulProp.plotVar = ['footprint']

# setting up mesh extension options
simulProp.meshExtensionAllDir = True
simulProp.maxElementIn = 10000
simulProp.set_mesh_extension_factor(1.1)
simulProp.set_mesh_extension_direction(['all'])
simulProp.meshReductionPossible = True
simulProp.maxCellSize = 0.1

# initialization parameters
Fr_geometry = Geometry('radial', radius=0.1, center=[0.1,0.])
init_param = InitializationParameters(Fr_geometry, regime='M')

# creating fracture object
Fr = Fracture(Mesh,
              init_param,
              Solid,
              Fluid,
              Injection,
              simulProp)

################################################################################
# the following lines are needed if you want to restart an existing simulation #
################################################################################
# from visualization import *
# Fr_list, properties = load_fractures(address="./Data/toughness_jump_positive", step_size=10)       # load all fractures                                                # list of times
# Solid, Fluid, Injection, simulProp = properties
# Fr = Fr_list[-1]

# create a Controller
controller = Controller(Fr,
                        Solid,
                        Fluid,
                        Injection,
                        simulProp)

# run the simulation
controller.run()

####################
# plotting results #
####################

from visualization import *

# loading simulation results
Fr_list, properties = load_fractures(address=myfolder,step_size=2)       # load all fractures
time_srs = get_fracture_variable(Fr_list, variable='time')                                                  # list of times
Solid, Fluid, Injection, simulProp = properties

# plot fracture radius
plot_prop = PlotProperties()
Fig_R = plot_fracture_list(Fr_list,
                           variable='footprint',
                           plot_prop=plot_prop)
Fig_R = plot_fracture_list(Fr_list,
                           fig=Fig_R,
                           variable='mesh',
                           mat_properties=properties[0],
                           backGround_param='K1c',
                           plot_prop=plot_prop)

# plot fracture radius
plot_prop = PlotProperties()
plot_prop.lineStyle = '.'               # setting the linestyle to point
plot_prop.graphScaling = 'loglog'       # setting to log log plot
Fig_R = plot_fracture_list(Fr_list,
                           variable='d_mean',
                           plot_prop=plot_prop)

# plot analytical radius
Fig_R = plot_analytical_solution(regime='M',
                                 variable='d_mean',
                                 mat_prop=Solid,
                                 inj_prop=Injection,
                                 fluid_prop=Fluid,
                                 time_srs=time_srs,
                                 fig=Fig_R)

Fig_R = plot_analytical_solution(regime='K',
                                 variable='d_mean',
                                 mat_prop=Solid,
                                 inj_prop=Injection,
                                 fluid_prop=Fluid,
                                 time_srs=time_srs,
                                 fig=Fig_R)


#  set block=True and comment last 2 lines if you want to keep the window open
#plt.show(block=True)
# plt.show(block=False)
# plt.pause(5)
# plt.close()


####################################
# exporting to multiple json files #
####################################

export_results = False

if export_results:


    # 0) import some functions needed later
    from visualization import *
    from postprocess_fracture import append_to_json_file


    # 1) decide the names of the Json files:
    myJsonName_1 = "./Data/toughness_heterog_3p6_export.json"           # I will export here most of the infos
    myJsonName_2 = "./Data/toughness_heterog_3p6_VEL_as_vector.json"       # I will export here the velocity infos



    # 2) load the results:
    #
    # >>> Remember that you can select a subset of time steps <<<<
    # >>> otherwise you will export at all the time steps     <<<<
    #
    print("\n 1) loading results")
    Fr_list, properties = load_fractures(address=myfolder) # or load_fractures(address=myfolder,time_srs=np.linspace(5., 8.0,600))
    Solid, Fluid, Injection, simulProp = properties
    print(" <-- DONE\n")


    # *) write to json the general informations
    simul_info = {'Eprime': Solid.Eprime,
                  'max_KIc': Solid.K1c.max(),
                  'min_KIc': Solid.K1c.min(),
                  'max_Sigma0': Solid.SigmaO.max(),
                  'min_Sigma0': Solid.SigmaO.min(),
                  'viscosity': Fluid.viscosity,
                  'total_injection_rate': Injection.injectionRate.max(),
                  'sources_coordinates_lastFR': Fr_list[-1].mesh.CenterCoor[Injection.sourceElem].tolist()}
    append_to_json_file(myJsonName_1, simul_info, 'append2keyASnewlist', key='simul_info',delete_existing_filename=True) #be careful: delete_existing_filename=True only the first time you call "append_to_json_file"


    # 3) write to json the coordinates of the points defining the fracture fronts at each time:
    print("\n 2) writing fronts")
    time_srs = get_fracture_variable(Fr_list,variable='time') # get the list of times corresponding to each fracture object
    append_to_json_file(myJsonName_1, time_srs, 'append2keyASnewlist', key='time_srs_of_Fr_list')
    fracture_fronts = []
    numberof_fronts = [] #there might be multiple fracture fronts in general
    mesh_info = [] # if you do not make remeshing or mesh extension you can export it only once
    index = 0
    for fracture in Fr_list:
        fracture_fronts.append(np.ndarray.tolist(fracture.Ffront))
        numberof_fronts.append(fracture.number_of_fronts)
        mesh_info.append([Fr_list[index].mesh.Lx, Fr_list[index].mesh.Ly, Fr_list[index].mesh.nx, Fr_list[index].mesh.ny])
        index = index + 1
    append_to_json_file(myJsonName_1, fracture_fronts, 'append2keyASnewlist', key='Fr_list')
    append_to_json_file(myJsonName_1, numberof_fronts, 'append2keyASnewlist', key='Number_of_fronts')
    append_to_json_file(myJsonName_1,mesh_info,'append2keyASnewlist', key='mesh_info')
    print(" <-- DONE\n")



    # 4) get the fracture opening w(t) versus time (t) at a point of coordinates myX and myY
    print("\n 3) get w(t) at a point... ")
    my_X = 0.02 ; my_Y = 0.
    w_at_my_point, time_list_at_my_point = get_fracture_variable_at_point(Fr_list, variable='w', point=[my_X, my_Y])
    append_to_json_file(myJsonName_1, w_at_my_point, 'append2keyASnewlist', key='w_at_my_point')
    append_to_json_file(myJsonName_1, time_list_at_my_point, 'append2keyASnewlist', key='time_list_W_at_my_point')
    print(" <-- DONE\n")



    # 5) get the fluid pressure p(t) versus time (t) at a point of coordinates myX and myY
    print("\n 4) get pf(t) at a point... ")
    my_X = 0.02 ; my_Y = 0.
    pf_at_my_point, time_list_at_my_point = get_fracture_variable_at_point(Fr_list, variable='pf', point=[my_X, my_Y])
    append_to_json_file(myJsonName_1, pf_at_my_point, 'append2keyASnewlist', key='pf_at_my_point_A')
    append_to_json_file(myJsonName_1, time_list_at_my_point, 'append2keyASnewlist', key='time_list_pf_at_my_point_A')
    print(" <-- DONE\n")


    # 5) get the fluid pressure p(t) versus time (t) at a point of coordinates myX and myY
    print("\n 4) get pf(t) at a point... ")
    my_X = 0.00 ; my_Y = 0.00
    pf_at_my_point, time_list_at_my_point = get_fracture_variable_at_point(Fr_list, variable='pf', point=[my_X, my_Y])
    append_to_json_file(myJsonName_1, pf_at_my_point, 'append2keyASnewlist', key='pf_at_my_point_B')
    append_to_json_file(myJsonName_1, time_list_at_my_point, 'append2keyASnewlist', key='time_list_pf_at_my_point_B')
    print(" <-- DONE\n")


    # 6) get w(y) along a vertical line passing through mypoint for different times
    print("\n 5) get w(y) with y passing through a specific point for different times... ")
    my_X = 0.; my_Y = 0.
    ext_pnts = np.empty((2, 2), dtype=np.float64)
    fracture_list_slice = plot_fracture_list_slice(Fr_list,
                                                   variable='w',
                                                   projection='2D',
                                                   plot_cell_center=True,
                                                   extreme_points=ext_pnts,
                                                   orientation='horizontal',
                                                   point1=[my_X , my_Y],
                                                   export2Json=True,
                                                   export2Json_assuming_no_remeshing=False)
    towrite = {'w_vert_slice_': fracture_list_slice}
    append_to_json_file(myJsonName_1, towrite, 'extend_dictionary')
    print(" <-- DONE\n")



    # 7) get pf(x) along a horizontal line passing through mypoint for different times
    print("\n 6) get pf(x) with x passing through a specific point for different times... ")
    my_X = 0.; my_Y = 0.
    ext_pnts = np.empty((2, 2), dtype=np.float64)
    fracture_list_slice = plot_fracture_list_slice(Fr_list,
                                                   variable='pf',
                                                   projection='2D',
                                                   plot_cell_center=True,
                                                   extreme_points=ext_pnts,
                                                   orientation='horizontal',
                                                   point1=[my_X , my_Y],
                                                   export2Json=True,
                                                   export2Json_assuming_no_remeshing=False)
    towrite = {'pf_horiz_slice_': fracture_list_slice}
    append_to_json_file(myJsonName_1, towrite, 'extend_dictionary')
    print(" <-- DONE\n")



    # 8) get ux(x) along a horizontal line passing through mypoint for different times
    print("\n 7) writing ux, qx")
    from postprocess_fracture import get_velocity_slice
    my_X = 0.01; my_Y = 0.
    ux_val, ux_times, ux_coord_x = get_velocity_slice(Solid, Fluid, Fr_list, [my_X, my_Y], orientation='horizontal')
    towrite = {'ux_horizontal_y0_value': ux_val,
               'ux_horizontal_y0_time': ux_times,
               'ux_horizontal_y0_coord': ux_coord_x}
    append_to_json_file(myJsonName_1, towrite, 'extend_dictionary')
    print(" <-- DONE\n")



    # 9) get uy(y) along a vertical line passing through mypoint for different times
    print("\n 8) writing uy, qy")
    from postprocess_fracture import get_velocity_slice
    my_X = 0.01; my_Y = 0.
    uy_val, uy_times, uy_coord_y = get_velocity_slice(Solid, Fluid, Fr_list, [my_X, my_Y], vel_direction='uy',
                                                      orientation='vertical')
    towrite = {'uy_vertical_x0_value': uy_val,
               'uy_vertical_x0_time': uy_times,
               'uy_vertical_x0_coord': uy_coord_y}
    append_to_json_file(myJsonName_1, towrite, 'extend_dictionary')
    print(" <-- DONE\n")



    # 10) get w(x,y,t) and pf(x,y,t)
    print("\n 9) get w(x,y,t) and  pf(x,y,t)... ")
    wofxyandt = []
    pofxyandt = []
    info = []
    jump = True #this is used to jump the first fracture
    for frac in Fr_list:
        if not jump:
            wofxyandt.append(np.ndarray.tolist(frac.w))
            pofxyandt.append(np.ndarray.tolist(frac.pFluid))
            info.append([frac.mesh.Lx,frac.mesh.Ly,frac.mesh.nx,frac.mesh.ny,frac.time])
        else:
            jump = False

    append_to_json_file(myJsonName_1, wofxyandt, 'append2keyASnewlist', key='w')
    append_to_json_file(myJsonName_1, pofxyandt, 'append2keyASnewlist', key='p')
    append_to_json_file(myJsonName_1, info, 'append2keyASnewlist', key='info_for_w_and_p')
    print(" <-- DONE\n")



    # 10) get fields of velocity and fluxes at selected times
    print("\n 10) process v(x,y), q(x,y)")
        # NOTE: saving only the non-zero-entries
    selected_times = range(len(Fr_list))
    vel_times = []   # list of the selected times
    vel_fields = []  # list with non zero velocity components (for each selected time)
    flux_fields = [] # list with non zero flux components (for each selected time)
    mesh_info = []   # list with the mesh info for each time
    index = 0

    for fracture in Fr_list: # loop on the fractures and take the data only from the selected indexes
        if index != 0 and index in selected_times:
            localVlist = []
            for i in range(fracture.fluidVelocity_components.shape[1]):
                localElementList = np.ndarray.tolist(fracture.fluidVelocity_components[:, i])
                if fracture.fluidVelocity_components[:, i].max() != 0 and i in fracture.EltChannel:
                    localElementList.append(i)
                    localVlist.append(localElementList)
            vel_fields.append(localVlist)
            vel_times.append(fracture.time)
            mesh_info.append([Fr_list[index].mesh.Lx, Fr_list[index].mesh.Ly, Fr_list[index].mesh.nx, Fr_list[index].mesh.ny])
            # flux_fields.append(np.ndarray.tolist(fracture.fluidFlux_components))
        index = index + 1

    append_to_json_file(myJsonName_2, vel_fields, 'append2keyASnewlist', key='vel_list', delete_existing_filename=True)
    append_to_json_file(myJsonName_2, vel_times, 'append2keyASnewlist', key='vel_times')
    append_to_json_file(myJsonName_2, mesh_info, 'append2keyASnewlist', key='mesh_info')
    # append_to_json_file(myJsonName_1, flux_fields, 'append2keyASnewlist', key='flux_list')
    print(" <-- DONE\n")

    print("DONE! in " + myJsonName_1)
    print("DONE! in " + myJsonName_2)
    plt.show(block=True)