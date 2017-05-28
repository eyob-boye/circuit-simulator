import sys
from django.shortcuts import render
from django.http import HttpResponse
from django.forms.models import model_to_dict
from models import SimulationCase, SimulationCaseForm, CircuitSchematics, CircuitSchematicsForm, CircuitComponents
import models
import os
import network_reader as NwRdr
import circuit_elements as CktElem
import solver as Slv
import matrix
import multiprocessing
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as Matpyplot


def prepare_simulation_objects(sim_para_model, sim_components, ckt_topo, conn_ckt_mat):
    
    synthesized_ckt_comps = {}
    
    components_found = sim_components[0]
    component_objects = sim_components[1]

    synthesized_ckt_comps["components_found"] = components_found
    synthesized_ckt_comps["component_objects"] = component_objects

    node_list = ckt_topo[0]
    branch_map = ckt_topo[1]

    synthesized_ckt_comps["branch_map"] = branch_map
    synthesized_ckt_comps["node_list"] = node_list

    # Segregating the elements in the circuit as those
    # having source voltages, meters and those that can be
    # controlled.
    bundled_list_of_elements = NwRdr.classify_components(components_found, component_objects)
    
    source_list = bundled_list_of_elements[0]
    meter_list = bundled_list_of_elements[1]
    controlled_elements = bundled_list_of_elements[2]
    
    synthesized_ckt_comps["source_list"] = source_list
    synthesized_ckt_comps["meter_list"] = meter_list
    synthesized_ckt_comps["controlled_elements"] = controlled_elements

    # Make a list of all the loops in the circuit.
    loop_list, loop_branches = NwRdr.determine_loops(conn_ckt_mat, node_list, branch_map)

    synthesized_ckt_comps["loop_list"] = loop_list
    synthesized_ckt_comps["loop_branches"] = loop_branches

    # Convert the above list of loops as segments of branches
    # while branch_params lists out the branch segments with
    # their parameters as the last element
    system_loops, branch_params = NwRdr.update_branches_loops(loop_branches, source_list)

    synthesized_ckt_comps["system_loops"] = system_loops
    synthesized_ckt_comps["branch_params"] = branch_params

    # Make lists of all the nodes connected together by 
    # empty branches.
    shortnode_list, shortbranch_list = NwRdr.delete_empty_branches(system_loops, \
                    branch_params, node_list, component_objects)

    synthesized_ckt_comps["shortnode_list"] = shortnode_list
    synthesized_ckt_comps["shortbranch_list"] = shortbranch_list

    # Make a list which has the components of a branch
    # as every element of the list.
    components_in_branch = []
    for c1 in range(len(branch_params)):
        current_branch_vector = []
        for c2 in range(len(branch_params[c1][:-1])):
            try:
                comp_pos = NwRdr.csv_element(branch_params[c1][c2])
                component_objects[comp_pos]
            except:
                pass
            else:
                current_branch_vector.append(comp_pos)
        components_in_branch.append(current_branch_vector)

    synthesized_ckt_comps["components_in_branch"] = components_in_branch

    # Branch currents for nodal analysis
    branch_currents = []
    for c1 in range(len(branch_params)):
        branch_currents.append(0.0)

    synthesized_ckt_comps["branch_currents"] = branch_currents

    # This is a list of branches which are only
    # in stiff loops. This comes in handy to recalculate
    # currents through branches that have recently become
    # stiff.
    branches_in_stiff_loops = []
    for c1 in range(len(branch_params)):
        branches_in_stiff_loops.append("yes")

    synthesized_ckt_comps["branches_in_stiff_loops"] = branches_in_stiff_loops

    # Stiff ratio indicates whether a branch
    # is stiff.
    stiff_ratio = []
    for c1 in range(len(branch_params)):
        stiff_ratio.append("no")

    synthesized_ckt_comps["stiff_ratio"] = stiff_ratio

    # An event vector which corresponds to every branch
    # in the circuit. An event is generated if any element
    # in the branch changes. Intitally, an event is generated
    # for each branch in the circuit.
    branch_events = []
    for c1 in range(len(branch_params)):
        branch_events.append("yes")

    synthesized_ckt_comps["branch_events"] = branch_events

    # Saving the previous set of branch events
    # to look for device state changes.
    branch_events_prev = []
    for c1 in range(len(branch_params)):
        branch_events_prev.append("no")

    synthesized_ckt_comps["branch_events_prev"] = branch_events_prev

    # Mark which loops are stiff.
    loop_stiff_info = []
    for c1 in range(len(system_loops)):
        loop_stiff_info.append("no")

    synthesized_ckt_comps["loop_stiff_info"] = loop_stiff_info

    # Node voltages for nodal analysis
    node_voltage = []
    for c1 in range(len(node_list)):
        node_voltage.append(0.0)

    synthesized_ckt_comps["node_voltage"] = node_voltage

    # Collect branches into bundles - those with
    # inductances, those with nonlinear elements,
    # those with voltage sources.
    bundled_list_of_branches = NwRdr.classify_branches(branch_params, component_objects)

    nonlinear_freewheel_branches = bundled_list_of_branches[0]
    inductor_list = bundled_list_of_branches[1][0]
    inductor_stiffness = bundled_list_of_branches[1][1]
    voltmeter_branches = bundled_list_of_branches[2][0]
    voltmeter_voltages = bundled_list_of_branches[2][1]

    synthesized_ckt_comps["nonlinear_freewheel_branches"] = nonlinear_freewheel_branches
    synthesized_ckt_comps["inductor_list"] = inductor_list
    synthesized_ckt_comps["inductor_stiffness"] = inductor_stiffness
    synthesized_ckt_comps["voltmeter_branches"] = voltmeter_branches
    synthesized_ckt_comps["voltmeter_voltages"] = voltmeter_voltages

    # For debugging, lists out a branch as a collection of 
    # the elements in the branch.
    branch_tags_in_loops = NwRdr.human_branch_names(branch_params, component_objects)

    synthesized_ckt_comps["branch_tags_in_loops"] = branch_tags_in_loops

    # These are lists needed to perform the KCL. Contains 
    # all the branches connected at a branch
    kcl_node_list,  abridged_node_voltage,  kcl_branch_map,  \
            branches_in_kcl_nodes,  admittance_matrix,  source_vector = \
            NwRdr.determine_kcl_parameters(branch_params,  node_list,  shortnode_list)

    synthesized_ckt_comps["kcl_node_list"] = kcl_node_list
    synthesized_ckt_comps["abridged_node_voltage"] = abridged_node_voltage
    synthesized_ckt_comps["kcl_branch_map"] = kcl_branch_map
    synthesized_ckt_comps["branches_in_kcl_nodes"] = branches_in_kcl_nodes
    synthesized_ckt_comps["admittance_matrix"] = admittance_matrix
    synthesized_ckt_comps["source_vector"] = source_vector
    
    # These are dictionaries to take snapshops of the system
    # at every event. The system loop map, and loop current
    # calculation information is stored when a new event is
    # detected and is retrieved when the event repeats.
    snap_branch_stiffness = {}
    snap_system_loopmap = {}
    snap_nonstiff_loops = {}
    snap_single_collection_nonstiff = {}
    snap_compute_loops_nonstiff = {}
    snap_loop_map_collection_nonstiff = {}
    snap_single_collection_stiff = {}
    snap_compute_loops_stiff = {}
    snap_loop_map_collection_stiff = {}

    synthesized_ckt_comps["snap_branch_stiffness"] = snap_branch_stiffness
    synthesized_ckt_comps["snap_system_loopmap"] = snap_system_loopmap
    synthesized_ckt_comps["snap_nonstiff_loops"] = snap_nonstiff_loops
    synthesized_ckt_comps["snap_single_collection_nonstiff"] = snap_single_collection_nonstiff
    synthesized_ckt_comps["snap_compute_loops_nonstiff"] = snap_compute_loops_nonstiff
    synthesized_ckt_comps["snap_loop_map_collection_nonstiff"] = snap_loop_map_collection_nonstiff
    synthesized_ckt_comps["snap_single_collection_stiff"] = snap_single_collection_stiff
    synthesized_ckt_comps["snap_compute_loops_stiff"] = snap_compute_loops_stiff
    synthesized_ckt_comps["snap_loop_map_collection_stiff"] = snap_loop_map_collection_stiff

    # Initialize the system loop map.
    system_loop_map = []
    synthesized_ckt_comps["system_loop_map"] = system_loop_map

    # System matrices for the ODE
    system_size = len(loop_branches)
    synthesized_ckt_comps["system_size"] = system_size

    sys_mat_a = matrix.Matrix(system_size, system_size)
    sys_mat_e = matrix.Matrix(system_size, system_size)
    curr_state_vec = matrix.Matrix(system_size)
    next_state_vec = matrix.Matrix(system_size)

    synthesized_ckt_comps["sys_mat_a"] = sys_mat_a
    synthesized_ckt_comps["sys_mat_e"] = sys_mat_e
    synthesized_ckt_comps["curr_state_vec"] = curr_state_vec
    synthesized_ckt_comps["next_state_vec"] = next_state_vec

    # These vectors are for the reduced order KVL
    reduced_curr_state = matrix.Matrix(system_size)
    reduced_next_state = matrix.Matrix(system_size)

    synthesized_ckt_comps["reduced_curr_state"] = reduced_curr_state
    synthesized_ckt_comps["reduced_next_state"] = reduced_next_state

    if source_list:
        source_size = len(source_list)
        sys_mat_b = matrix.Matrix(system_size, len(source_list))
        sys_mat_u = matrix.Matrix(len(source_list))

    else:
        source_size = 1
        sys_mat_b = matrix.Matrix(system_size)
        sys_mat_u = 0.0

    synthesized_ckt_comps["source_size"] = source_size
    synthesized_ckt_comps["sys_mat_b"] = sys_mat_b
    synthesized_ckt_comps["sys_mat_u"] = sys_mat_u

    # 5th order Runge Kutta method
    ##ode_k1=matrix.Matrix(system_size)
    ##ode_k2=matrix.Matrix(system_size)
    ##ode_k3=matrix.Matrix(system_size)
    ##ode_k4=matrix.Matrix(system_size)
    ##ode_k5=matrix.Matrix(system_size)
    ##ode_k6=matrix.Matrix(system_size)
    ##ode_dbydt=matrix.Matrix(system_size)
    ##ode_var=[ode_k1, ode_k2, ode_k3, ode_k4, ode_k5, ode_k6, ode_dbydt]

    # 4th order Runge Kutta method
    ode_k1 = matrix.Matrix(system_size)
    ode_k2 = matrix.Matrix(system_size)
    ode_k3 = matrix.Matrix(system_size)
    ode_k4 = matrix.Matrix(system_size)
    ode_dbydt = matrix.Matrix(system_size)
    ode_var = [ode_k1, ode_k2, ode_k3, ode_k4, ode_dbydt]

    synthesized_ckt_comps["ode_k1"] = ode_k1
    synthesized_ckt_comps["ode_k2"] = ode_k2
    synthesized_ckt_comps["ode_k3"] = ode_k3
    synthesized_ckt_comps["ode_k4"] = ode_k4
#    synthesized_ckt_comps["ode_k5"] = ode_k5
#    synthesized_ckt_comps["ode_k6"] = ode_k6
    synthesized_ckt_comps["ode_dbydt"] = ode_dbydt
    synthesized_ckt_comps["ode_var"] = ode_var

    # Trapezoidal rule
    ##ode_k1=matrix.Matrix(system_size)
    ##ode_k2=matrix.Matrix(system_size)
    ##ode_dbydt=matrix.Matrix(system_size)
    ##ode_var=[ode_k1, ode_k2, ode_dbydt]

    # Assign the parameters to the circuit elements.
#    NwRdr.read_circuit_parameters(nw_input, conn_ckt_mat, branch_params, component_objects, components_found)

#    # Read the control code and the descriptions of the control code
#    control_files,  control_functions,  control_file_inputs, \
#            control_file_outputs, control_file_staticvars, control_file_timeevents, \
#            control_file_variablestorage,  control_file_events = \
#            NwRdr.update_control_code(component_objects, components_found, \
#                                      controlled_elements, meter_list, control_files)

    control_files, control_functions, control_file_inputs, \
            control_file_outputs, control_file_staticvars, control_file_timeevents, \
            control_file_variablestorage, control_file_events = \
            NwRdr.update_control_django(sim_para_model, component_objects, components_found, \
                                      controlled_elements, meter_list)

    synthesized_ckt_comps["control_files"] = control_files
    synthesized_ckt_comps["control_functions"] = control_functions
    synthesized_ckt_comps["control_file_inputs"] = control_file_inputs
    synthesized_ckt_comps["control_file_outputs"] = control_file_outputs
    synthesized_ckt_comps["control_file_staticvars"] = control_file_staticvars
    synthesized_ckt_comps["control_file_timeevents"] = control_file_timeevents
    synthesized_ckt_comps["control_file_variablestorage"] = control_file_variablestorage
    synthesized_ckt_comps["control_file_events"] = control_file_events

    plotted_variable_list = []
    for c1 in control_file_variablestorage.keys():
        if control_file_variablestorage[c1][1].lower()=="yes" or \
                control_file_variablestorage[c1][1].lower()=="y":
            plotted_variable_list.append(c1)

    print(control_file_variablestorage)
    print
    print
    print(plotted_variable_list)
    print
    print
    synthesized_ckt_comps["plotted_variable_list"] = plotted_variable_list

    return synthesized_ckt_comps



def simulation_iterations(sim_id, synthesized_ckt_comps):
    sim_para_model = SimulationCase.objects.get(id=int(sim_id))
    t = 0.0
    t_ode = 0.0
    dt = sim_para_model.sim_time_step
    t_ode_prev = -sim_para_model.sim_time_step
    t_diff = sim_para_model.sim_time_step
    dt_store = sim_para_model.sim_time_data
    t_store = dt_store
    t_limit = sim_para_model.sim_time_limit
    
    outputfile = []
    outputfilename = sim_para_model.sim_output_file.split(".")[0]
    outputfileext = sim_para_model.sim_output_file.split(".")[1]
    if sim_para_model.sim_output_slice=="No":
        file_path = os.path.join(os.sep, \
                        sim_para_model.sim_working_directory, \
                        sim_para_model.sim_output_file)
        outputfile.append(open(file_path, "w"))
    else:
        for c1 in range(sim_para_model.sim_div_number):
            file_path = os.path.join(os.sep, \
                                        sim_para_model.sim_working_directory, \
                                        outputfilename+str(c1+1)+"."+outputfileext)
                                        
            file_object = open(file_path, "w")
            outputfile.append(file_object)
    
    t_split_file = []
    t_window = t_limit/sim_para_model.sim_div_number
    t_window_start = 0.0
    for c1 in range(sim_para_model.sim_div_number):
        t_split_file.append(t_window_start)
        t_window_start += t_window
    t_split_file.append(t_window_start)

    sys.path.insert(0, sim_para_model.sim_working_directory)
    import __control

    branch_events = synthesized_ckt_comps["branch_events"]
    branch_params = synthesized_ckt_comps["branch_params"]
    component_objects = synthesized_ckt_comps["component_objects"]
    source_list = synthesized_ckt_comps["source_list"]
    sys_mat_u = synthesized_ckt_comps["sys_mat_u"]
    components_in_branch = synthesized_ckt_comps["components_in_branch"]
    system_loops = synthesized_ckt_comps["system_loops"]
    system_loop_map = synthesized_ckt_comps["system_loop_map"]
    stiff_ratio = synthesized_ckt_comps["stiff_ratio"]

    branches_in_kcl_nodes = synthesized_ckt_comps["branches_in_kcl_nodes"]
    kcl_branch_map = synthesized_ckt_comps["kcl_branch_map"]
    branch_tags_in_loops = synthesized_ckt_comps["branch_tags_in_loops"]

    branch_events_prev = synthesized_ckt_comps["branch_events_prev"]
    snap_branch_stiffness = synthesized_ckt_comps["snap_branch_stiffness"]

    curr_state_vec = synthesized_ckt_comps["curr_state_vec"]
    next_state_vec = synthesized_ckt_comps["next_state_vec"]
    loop_stiff_info = synthesized_ckt_comps["loop_stiff_info"]
    snap_system_loopmap = synthesized_ckt_comps["snap_system_loopmap"]
    snap_single_collection_nonstiff = synthesized_ckt_comps["snap_single_collection_nonstiff"]
    snap_compute_loops_nonstiff = synthesized_ckt_comps["snap_compute_loops_nonstiff"]
    snap_loop_map_collection_nonstiff = synthesized_ckt_comps["snap_loop_map_collection_nonstiff"]
    snap_single_collection_stiff = synthesized_ckt_comps["snap_single_collection_stiff"]
    snap_compute_loops_stiff = synthesized_ckt_comps["snap_compute_loops_stiff"]
    snap_loop_map_collection_stiff = synthesized_ckt_comps["snap_loop_map_collection_stiff"]
    snap_nonstiff_loops = synthesized_ckt_comps["snap_nonstiff_loops"]
    sys_mat_e = synthesized_ckt_comps["sys_mat_e"]
    sys_mat_a = synthesized_ckt_comps["sys_mat_a"]
    sys_mat_b = synthesized_ckt_comps["sys_mat_b"]
    branch_currents = synthesized_ckt_comps["branch_currents"]
    admittance_matrix = synthesized_ckt_comps["admittance_matrix"]
    source_vector = synthesized_ckt_comps["source_vector"]
    abridged_node_voltage = synthesized_ckt_comps["abridged_node_voltage"]
    node_voltage = synthesized_ckt_comps["node_voltage"]
    kcl_node_list = synthesized_ckt_comps["kcl_node_list"]
    node_list = synthesized_ckt_comps["node_list"]
    shortnode_list = synthesized_ckt_comps["shortnode_list"]
    nonlinear_freewheel_branches = synthesized_ckt_comps["nonlinear_freewheel_branches"]
    reduced_curr_state = synthesized_ckt_comps["reduced_curr_state"]
    reduced_next_state = synthesized_ckt_comps["reduced_next_state"]
    ode_var = synthesized_ckt_comps["ode_var"]
    system_size = synthesized_ckt_comps["system_size"]
    meter_list = synthesized_ckt_comps["meter_list"]

    plotted_variable_list = synthesized_ckt_comps["plotted_variable_list"]
    control_files = synthesized_ckt_comps["control_files"]
    control_files = synthesized_ckt_comps["control_files"]
    control_functions = synthesized_ckt_comps["control_functions"]
    control_file_inputs = synthesized_ckt_comps["control_file_inputs"]
    control_file_outputs = synthesized_ckt_comps["control_file_outputs"]
    control_file_staticvars = synthesized_ckt_comps["control_file_staticvars"]
    control_file_timeevents = synthesized_ckt_comps["control_file_timeevents"]
    control_file_variablestorage = synthesized_ckt_comps["control_file_variablestorage"]
    control_file_events = synthesized_ckt_comps["control_file_events"]
#    control_files, control_functions, control_file_inputs, \
#            control_file_outputs, control_file_staticvars, \
#            control_file_timeevents, control_file_variablestorage, control_file_events

#    control_file_inputs, control_file_outputs, control_file_staticvars, control_file_timeevents, \
#                            control_file_variablestorage, control_file_events, component_objects


    while t<t_limit:

        # Check if an event has been generated. If yes,
        # recalculate the system matrices. If not, just
        # solve the ODE
        
        if ("yes" in branch_events) or ("hard" in branch_events):
            
            # Initialize the branch parameters from the
            # values of the components in the branch.

            NwRdr.initialize_branch_params(branch_params, branch_events, component_objects, source_list, sys_mat_u, components_in_branch)
            
            # This system map is to describe whether in each loop
            # a branch exists - if it does, is it "stiff" or not ("yes")
            # Initialize all of the elements to "no"
            if not system_loop_map:
                for c1 in range(len(system_loops)):
                    br_vector = []
                    for c2 in range(len(branch_params)):
                        br_vector.append("no")
                    system_loop_map.append(br_vector)
    
                # This function will populate the system loop map from the system loops
                # Also, this will create the list stiff_ratio that defines which
                # branches are stiff
                Slv.generate_system_loops(system_loops, branch_params, system_loop_map, stiff_ratio, dt)

                # Delete any loops that violate the basic rules
                for c1 in range(len(system_loop_map)-1, -1, -1):
                    is_loop_genuine = Slv.loop_validity_checking(system_loop_map[c1], branches_in_kcl_nodes, kcl_branch_map)
                    if is_loop_genuine=="no":
                        del system_loop_map[c1]
                

                # The next block of code seeks to eliminate redundant loops
                # or dependent loops. These are loops that are not identical to
                # other loops but are linear combinations of other loops. If not
                # eliminated, these extra loops can cause instability.
                
                # The loops are arranged into clusters according to the first
                # branch in these loops.
                loop_clusters = []
                c1 = 0
                c2 = 0
                while (c1<len(system_loop_map)) and (c2<len(system_loop_map[0])):
                    # For each loop, attempt to loop for the first branch.
                    curr_loop_cluster = []
                    
                    # Assume that the first branch c2 in the loop c1
                    # is not a branch. In that case loop for subsequent
                    # loops to check if any of those loops have branches
                    # at c2.
                    branch_found = "no"
                    if system_loop_map[c1][c2]=="no":
                        c3 = c1 + 1
                        c4 = c2
                        while c3<len(system_loop_map) and c4<len(system_loop_map[c1]) and branch_found=="no":
                            if not system_loop_map[c3][c4]=="no":
                                # Interchange the loops if another loop below it
                                # has a branch.
                                branch_found = "yes"
                                if not c1==c3:
                                    system_loop_map[c1], system_loop_map[c3] = system_loop_map[c3], system_loop_map[c1]
                                c2 = c4
                            else:
                                c3 += 1
                                # If all the loops have been exhausted, it means the
                                # branch is only found in previous loops. So move on
                                # to the next branch and go back to the loop c1.
                                if c3==len(system_loop_map):
                                    c3 = c1
                                    c4 += 1
                    
                    else:
                        branch_found = "yes"
                    
                    # If a branch has been found, eliminate branches from 
                    # all subsequent loops. Sometimes, this is not possible
                    # as the loops are not compatible. In that case the cluster
                    # is expanded and loops are interchanged.
                    if branch_found=="yes":
                        curr_loop_cluster.append(c1)
                        # Start with the first loop after the cluster.
                        c3 = curr_loop_cluster[-1] + 1
                        while c3<len(system_loop_map):
                            # Attempt to perform manipulations with every 
                            # loop in the cluster. Once successful, exit
                            c4 = 0
                            loop_manip_success = "no"
                            common_branch_found = "no"
                            while c4<len(curr_loop_cluster) and common_branch_found=="no":
                                row_position = curr_loop_cluster[c4]
                                if not system_loop_map[c3][c2]=="no":
                                    if not system_loop_map[row_position][c2]=="no":
                                        common_branch_found = "yes"
                                        if system_loop_map[c3][c2]==system_loop_map[row_position][c2]:
                                            loop_manip_result = Slv.loop_manipulations(system_loop_map, branches_in_kcl_nodes, kcl_branch_map, \
                                                    row_position, c3, "difference", branch_params, branch_tags_in_loops)
                                        else:
                                            loop_manip_result = Slv.loop_manipulations(system_loop_map, branches_in_kcl_nodes, kcl_branch_map, \
                                                    row_position, c3, "addition", branch_params, branch_tags_in_loops)
                                        
                                        if loop_manip_result:
                                            loop_manip_success = "yes"
                                
                                # Look for a common branch with every loop in the cluster
                                if common_branch_found=="no":
                                    c4 += 1
                                
                            # If a common branch has been found, but no loop
                            # manipulation was successful, the loop will have to
                            # be added to the cluster. 
                            if loop_manip_success=="no" and common_branch_found=="yes":
                                # If the loop is away from the last loop in the cluster,
                                # a loop interchange needs to be performed to bring the loop
                                # just below the cluster
                                if c3>curr_loop_cluster[-1]+1:
                                    system_loop_map[c3], system_loop_map[curr_loop_cluster[-1]+1] = system_loop_map[curr_loop_cluster[-1]+1], \
                                                    system_loop_map[c3]
                                    # Increment the loop pointer as the new loop will come below
                                    # the cluster.
                                    c3 = curr_loop_cluster[-1]+2
                                    curr_loop_cluster.append(c3-1)
                                else:
                                    c3 += 1
                                    curr_loop_cluster.append(c3)
                            else:
                                c3 += 1
                        
                        # Add the custer to all the loop clusters
                        loop_clusters.append(curr_loop_cluster)
                    
                    # Increment the loop and branch counters.
                    # The loop counter will be loop subsequent to the
                    # last loop in the last cluster.
#                    if loop_clusters:
                    c1 = loop_clusters[-1][-1] + 1
                    c2 += 1


                # Determine which branch each component is in. This is to 
                # speed up the simulator so as it does not search for a 
                # component in the circuit every time.
                for c1 in range(len(branch_params)):
                    for c2 in range(len(branch_params[c1][:-1])):
                        try:
                            comp_pos = NwRdr.csv_element(branch_params[c1][c2])
                            component_objects[comp_pos]
                        except:
                            pass
                        else:
                            component_objects[comp_pos].determine_branch(branch_params)

                
#                # No branch in the circuit that has an element in it can 
#                # have zero resistance. An empty branch with only "wire"
#                # is possible.
#                for c1 in range(len(branch_params)):
#                    if (branch_params[c1][-1][0][1] or \
#                        ((1.0 in branch_params[c1][-1][1]) or \
#                        (-1.0 in branch_params[c1][-1][1]))):
#                            
#                            if not branch_params[c1][-1][0][0]:
#                                
#                                print("Following branch has zero resistance")
#                                NwRdr.human_loop(branch_params[c1][:-1])
#                                print
#                                raise CktEx.BranchZeroResistanceError

                
                # This is to speed up the KCL calculations. This way the 
                # branches connected to nodes can be looked up from this
                # table and don't have to be calculated every time. The 
                # only branches that will be calculated are the nonlinear
                # branches.
                kcl_branch_lookup = []
                for c1 in range(len(branch_params)):
                    branch_vector = []
                    branch_vector.append(branch_events[c1])
                    for c2 in range(len(kcl_branch_map)):
                        for c3 in range(c2+1, len(kcl_branch_map[c2])):
                            if kcl_branch_map[c2][c3]:
                                if c1 in kcl_branch_map[c2][c3][0]:
                                    branch_vector.append([c2, c3])
                                    br_pos = kcl_branch_map[c2][c3][0].index(c1)
                                    branch_vector.append(kcl_branch_map[c2][c3][1][br_pos])
                        
                    if stiff_ratio[c1]=="yes":
                        branch_vector.append(1)
                    else:
                        if (branch_params[c1][-1][0][1]):
                            if  abs(branch_params[c1][-1][0][1]/branch_params[c1][-1][0][0]) < 0.1*dt:
                                branch_vector.append(1)
                            else:
                                branch_vector.append(2)
                        else:
                            branch_vector.append(1)
                    
                    branch_vector.append(branch_params[c1][-1][3])
                    branch_vector.append(branch_params[c1][-1][0][0])
                    branch_vector.append(branch_params[c1][-1][2])
                    kcl_branch_lookup.append(branch_vector)

            
            # This function will update the system loop map with
            # the stiffness information from previous iteration.
            Slv.update_system_loop_map(branch_params, system_loop_map, stiff_ratio, dt)
            
            # Store the branch events
            for c1 in range(len(branch_params)):
                branch_events_prev[c1] = branch_events[c1]


            # Each case of the simulation with a different set of
            # branch stiffness, will be stored for future reference to
            # reduce repeated calculations. These cases will be stored
            # in dictionaries and the keys will be the stiffness
            # information of the branches.
            branch_stiffness_key = ""
            for c1 in range(len(stiff_ratio)):
                if stiff_ratio[c1]=="yes":
                    branch_stiffness_key += "y"
                else:
                    branch_stiffness_key += "n"
            
            if branch_stiffness_key in snap_branch_stiffness.keys():
                branch_layout_found = "yes"
            else:
                branch_layout_found = "no"

            # Generate the system information matrices either by calculation or
            # by extracting the stored information from snapshots.
            
            if branch_layout_found=="no":
                # If the snapshot of the system interms of branch stiffness has not been
                # found, it means, that this state of the system has not been encountered
                # before. So all snapshot information is stored.
                
                # Check if the system is stiff and recalculate if it is.
                # The stiff loops are isolated to the minimum number of loops
                # So that the system dynamics are captured as far as possible.
                Slv.remove_stiffness(system_loop_map, [curr_state_vec, next_state_vec], loop_stiff_info, branches_in_kcl_nodes, kcl_branch_map, \
                                        branch_params, branch_tags_in_loops)

                Slv.approximate_nonstiff_loops(branch_params, stiff_ratio, system_loop_map, branches_in_kcl_nodes, kcl_branch_map)
                
                
                # With the stiff branches minimized to the minimum number
                # of loops, the stiffness information of the branches and the
                # manipulated system_loop_map are stored.
                snap_branch_stiffness[branch_stiffness_key] = []
                for c1 in range(len(stiff_ratio)):
                    snap_branch_stiffness[branch_stiffness_key].append(stiff_ratio[c1])
                
                snap_system_loopmap[branch_stiffness_key] = []
                for c1 in range(len(system_loop_map)):
                    row_vector = []
                    for c2 in range(len(system_loop_map[c1])):
                        row_vector.append(system_loop_map[c1][c2])
                    snap_system_loopmap[branch_stiffness_key].append(row_vector)
                

                single_nonstiff_collection = []
                compute_loops_nonstiff = []
                loop_map_collection_nonstiff = []
                single_stiff_collection = []
                compute_loops_stiff = []
                loop_map_collection_stiff = []
                nonstiff_loops = []
                
                
                
                # The next step is to divide the loops into stiff and non-stiff loops
                # After which information on how to compute the loop currents of the
                # non-stiff loops is stored in the snapshot dictionaries.
                single_nonstiff_collection, compute_loops_nonstiff, loop_map_collection_nonstiff, \
                single_stiff_collection, compute_loops_stiff, loop_map_collection_stiff, nonstiff_loops = \
                    Slv.compute_loop_currents_generate(branch_params, stiff_ratio, system_loop_map, branch_events, branches_in_kcl_nodes)

                
#                Slv.compute_loop_currents_calc(single_nonstiff_collection, compute_loops_nonstiff, loop_map_collection_nonstiff, \
#                            single_stiff_collection, compute_loops_stiff, loop_map_collection_stiff, branch_params, \
#                            [curr_state_vec, next_state_vec])
                
                # These are the branches that occur only once in a nonstiff loop
                # and so the branch currents automatically become the loop currents.
                snap_single_collection_nonstiff[branch_stiffness_key] = []
                for c1 in range(len(single_nonstiff_collection)):
                    row_vector = []
                    for c2 in range(len(single_nonstiff_collection[c1])):
                        row_vector.append(single_nonstiff_collection[c1][c2])
                    snap_single_collection_nonstiff[branch_stiffness_key].append(row_vector)

                # These are the nonstiff loops whose current needs to be calculated
                # as they do not have any unique branches.
                snap_compute_loops_nonstiff[branch_stiffness_key] = []
                for c1 in range(len(compute_loops_nonstiff)):
                    snap_compute_loops_nonstiff[branch_stiffness_key].append(compute_loops_nonstiff[c1])
                
                # The calculation of the loop currents is by solving equations AX=B
                # Below is the matrix A. The matrix B will always change as it is the
                # input.
                snap_loop_map_collection_nonstiff[branch_stiffness_key] = []
                for c1 in range(len(loop_map_collection_nonstiff)):
                    row_vector = []
                    for c2 in range(len(loop_map_collection_nonstiff[c1])):
                        row_vector.append(loop_map_collection_nonstiff[c1][c2])
                    snap_loop_map_collection_nonstiff[branch_stiffness_key].append(row_vector)

                # These are the branches that occur only once in a stiff loop
                # and so the branch currents automatically become the loop currents.
                snap_single_collection_stiff[branch_stiffness_key] = []
                for c1 in range(len(single_stiff_collection)):
                    row_vector = []
                    for c2 in range(len(single_stiff_collection[c1])):
                        row_vector.append(single_stiff_collection[c1][c2])
                    snap_single_collection_stiff[branch_stiffness_key].append(row_vector)

                # These are the stiff loops whose current needs to be calculated
                # as they do not have any unique branches.
                snap_compute_loops_stiff[branch_stiffness_key] = []
                for c1 in range(len(compute_loops_stiff)):
                    snap_compute_loops_stiff[branch_stiffness_key].append(compute_loops_stiff[c1])
                
                # The calculation of the loop currents is by solving equations AX=B
                # Below is the matrix A. The matrix B will always change as it is the
                # input.
                snap_loop_map_collection_stiff[branch_stiffness_key] = []
                for c1 in range(len(loop_map_collection_stiff)):
                    row_vector = []
                    for c2 in range(len(loop_map_collection_stiff[c1])):
                        row_vector.append(loop_map_collection_stiff[c1][c2])
                    snap_loop_map_collection_stiff[branch_stiffness_key].append(row_vector)
                
                # All the nonstiff loops. This is used for solving the final ODE.
                snap_nonstiff_loops[branch_stiffness_key] = []
                for c1 in range(len(nonstiff_loops)):
                    snap_nonstiff_loops[branch_stiffness_key].append(nonstiff_loops[c1])
                
            else:
                # If the system snapshot has been found before, all the information can be
                # immediately extracted from the snapshot matrices with minimal or no
                # calculations.
                
                system_loop_map = []
                for c1 in range(len(snap_system_loopmap[branch_stiffness_key])):
                    row_vector = []
                    for c2 in range(len(snap_system_loopmap[branch_stiffness_key][c1])):
                        row_vector.append(snap_system_loopmap[branch_stiffness_key][c1][c2])
                    system_loop_map.append(row_vector)
                

                single_nonstiff_collection = []
                for c1 in range(len(snap_single_collection_nonstiff[branch_stiffness_key])):
                    row_vector = []
                    for c2 in range(len(snap_single_collection_nonstiff[branch_stiffness_key][c1])):
                        row_vector.append(snap_single_collection_nonstiff[branch_stiffness_key][c1][c2])
                    single_nonstiff_collection.append(row_vector)
                

                compute_loops_nonstiff = []
                for c1 in range(len(snap_compute_loops_nonstiff[branch_stiffness_key])):
                    compute_loops_nonstiff.append(snap_compute_loops_nonstiff[branch_stiffness_key][c1])
                
                
                loop_map_collection_nonstiff = []
                for c1 in range(len(snap_loop_map_collection_nonstiff[branch_stiffness_key])):
                    row_vector = []
                    for c2 in range(len(snap_loop_map_collection_nonstiff[branch_stiffness_key][c1])):
                        row_vector.append(snap_loop_map_collection_nonstiff[branch_stiffness_key][c1][c2])
                    loop_map_collection_nonstiff.append(row_vector)

                single_stiff_collection = []
                for c1 in range(len(snap_single_collection_stiff[branch_stiffness_key])):
                    row_vector = []
                    for c2 in range(len(snap_single_collection_stiff[branch_stiffness_key][c1])):
                        row_vector.append(snap_single_collection_stiff[branch_stiffness_key][c1][c2])
                    single_stiff_collection.append(row_vector)
                

                compute_loops_stiff = []
                for c1 in range(len(snap_compute_loops_stiff[branch_stiffness_key])):
                    compute_loops_stiff.append(snap_compute_loops_stiff[branch_stiffness_key][c1])
                
                
                loop_map_collection_stiff = []
                for c1 in range(len(snap_loop_map_collection_stiff[branch_stiffness_key])):
                    row_vector = []
                    for c2 in range(len(snap_loop_map_collection_stiff[branch_stiffness_key][c1])):
                        row_vector.append(snap_loop_map_collection_stiff[branch_stiffness_key][c1][c2])
                    loop_map_collection_stiff.append(row_vector)

                nonstiff_loops = []
                for c1 in range(len(snap_nonstiff_loops[branch_stiffness_key])):
                    nonstiff_loops.append(snap_nonstiff_loops[branch_stiffness_key][c1])
            
            # When a branch becomes stiff naturally such as a diode or an IGBT turning off
            # because current through it reverses, it actual practice, the current is
            # extremely low and the reverse current that restores the charge is also extremely low.
            # In simulations however, when a device turns off the current can be significant.
            # When a branch becomes stiff for reason such as this, it should not result in
            # freewheeling as current is almost zero. Therefore, to ensure that freewheeling does
            # not take place, these branch currents are made zero.
            Slv.new_stiff_branch_adjustment(system_loop_map, branch_params, branch_events, stiff_ratio, [curr_state_vec, next_state_vec], \
                            sys_mat_e, sys_mat_a, sys_mat_b, sys_mat_u,  dt)

            
            # Calculate the inductor voltages
            # For debugging only
            #inductor_voltages = NwRdr.inductor_volt_calc(inductor_list, system_loop_map, branch_params, ode_var, dt)
            
            freewheel_attempt_no = 0
            new_branch_events = "yes"
            # This list is no longer used.
            nodes_in_kcl_calc = []
            while new_branch_events=="yes":
                
                # Initialize currents in branches
                # This is used to determine whether any
                # inductor current change event occurs.
                for c1 in range(len(branch_currents)):
                    branch_currents[c1] = 0.0

                # Update the look up table for the KCL calculations.
                for c1 in range(len(branch_params)):
                    kcl_branch_lookup[c1][0] = branch_events[c1]
#                    if branch_events[c1]=="yes" or branch_events[c1]=="hard":
                    if branch_events[c1]=="yes" or branch_events[c1]=="hard" or branch_params[c1][-1][0][1]:
                        if stiff_ratio[c1]=="yes":
                            kcl_branch_lookup[c1][3] = 1
                        else:
                            if (branch_params[c1][-1][0][1]):
                                if  abs(branch_params[c1][-1][0][1]/branch_params[c1][-1][0][0]) < 0.1*dt:
                                    kcl_branch_lookup[c1][3] = 1
                                else:
                                    kcl_branch_lookup[c1][3] = 2
                            else:
                                kcl_branch_lookup[c1][3] = 1

                    kcl_branch_lookup[c1][4] = branch_params[c1][-1][3]
                    kcl_branch_lookup[c1][5] = branch_params[c1][-1][0][0]
                    kcl_branch_lookup[c1][6] = branch_params[c1][-1][2]

                # Nodal analysis
                Slv.current_continuity(kcl_branch_lookup, admittance_matrix, source_vector, abridged_node_voltage, node_voltage, \
                                kcl_node_list, node_list, shortnode_list, branch_params, nonlinear_freewheel_branches, branch_currents, \
                                "det_state", t, 1)
                
                
                # Determine the state of nonlinear devices
                # Using the currents from nodal analysis,
                # it will check if any of the devices will start
                # or stop conducting.
                
                if freewheel_attempt_no > 1:
                    for comps in component_objects.keys():
                        component_objects[comps].determine_state(branch_currents, branch_params, branch_events)
                else:
                    for comps in component_objects.keys():
                        component_objects[comps].pre_determine_state(branch_currents, branch_params, branch_events)


                # Check if there is a difference in the branch events as opposed
                # to the previous nodal analysis. Keep performing nodal analysis
                # until there are no new branch events
                new_branch_events = "no"
                for c1 in range(len(branch_params)):
                    if not branch_events[c1]==branch_events_prev[c1]:
                        new_branch_events = "yes"

                
                # For those branches that have experienced an event,
                # get the new values of branch parameters from the
                # components in the branch.                    
                NwRdr.initialize_branch_params(branch_params, branch_events, component_objects, source_list, sys_mat_u, components_in_branch)


                # Store the branch events
                for c1 in range(len(branch_params)):
                    branch_events_prev[c1] = branch_events[c1]
                
                freewheel_attempt_no += 1
                
            # This is the end of the block to determine freewheeling.


            # This function will update the system loop map with
            # the stiffness information from continuity check.
            Slv.update_system_loop_map(branch_params, system_loop_map, stiff_ratio, dt)

           
            # Arrange the nonstiff loops in ascending order of di/dt.
##            Slv.arrange_nonstiff_loops(sys_mat_e, sys_mat_a, sys_mat_b, sys_mat_u,  dt, branch_params, stiff_ratio, system_loop_map,\
##                    [curr_state_vec, next_state_vec],  branch_events)


            # To obtain the system information matrices, check if the stiffness information
            # has been encountered before. If so, extract from the snapshot dictionaries
            # or else run the remove_stiffness and compute_loop_currents functions.
            branch_stiffness_key = ""
            for c1 in range(len(stiff_ratio)):
                if stiff_ratio[c1]=="yes":
                    branch_stiffness_key += "y"
                else:
                    branch_stiffness_key += "n"
            
            if branch_stiffness_key in snap_branch_stiffness.keys():
                branch_layout_found = "yes"
            else:
                branch_layout_found = "no"
            
            
            if branch_layout_found=="no":
                
                # Check if the system is stiff and recalculate if it is.
                # The stiff loops are isolated to the minimum number of loops
                # So that the system dynamics are captured as far as possible.
                Slv.remove_stiffness(system_loop_map, [curr_state_vec, next_state_vec], loop_stiff_info, branches_in_kcl_nodes, kcl_branch_map, \
                                        branch_params, branch_tags_in_loops)
                    
                # Add all the system matrices to the dictionary
                # corresponding to the layout of the branch stiff
                # information.
                
                Slv.approximate_nonstiff_loops(branch_params, stiff_ratio, system_loop_map, branches_in_kcl_nodes, kcl_branch_map)
                
                snap_branch_stiffness[branch_stiffness_key] = []
                for c1 in range(len(stiff_ratio)):
                    snap_branch_stiffness[branch_stiffness_key].append(stiff_ratio[c1])
                
                
                snap_system_loopmap[branch_stiffness_key] = []
                for c1 in range(len(system_loop_map)):
                    row_vector = []
                    for c2 in range(len(system_loop_map[c1])):
                        row_vector.append(system_loop_map[c1][c2])
                    snap_system_loopmap[branch_stiffness_key].append(row_vector)
                

                single_nonstiff_collection = []
                compute_loops_nonstiff = []
                loop_map_collection_nonstiff = []
                single_stiff_collection = []
                compute_loops_stiff = []
                loop_map_collection_stiff = []
                nonstiff_loops = []
                
                
                # The generate function generates the matrices for 
                # the loops that result from a given layout of branches.
                single_nonstiff_collection, compute_loops_nonstiff, loop_map_collection_nonstiff, \
                single_stiff_collection, compute_loops_stiff, loop_map_collection_stiff, nonstiff_loops = \
                    Slv.compute_loop_currents_generate(branch_params, stiff_ratio, system_loop_map, branch_events, branches_in_kcl_nodes)
                
                
#                Slv.compute_loop_currents_calc(single_collection_nonstiff, compute_loops_nonstiff, loop_map_collection_nonstiff, \
#                            single_collection_stiff, compute_loops_stiff, loop_map_collection_stiff, branch_params, \
#                            [curr_state_vec, next_state_vec])
                

                snap_single_collection_nonstiff[branch_stiffness_key] = []
                for c1 in range(len(single_nonstiff_collection)):
                    row_vector = []
                    for c2 in range(len(single_nonstiff_collection[c1])):
                        row_vector.append(single_nonstiff_collection[c1][c2])
                    snap_single_collection_nonstiff[branch_stiffness_key].append(row_vector)
                
                
                snap_compute_loops_nonstiff[branch_stiffness_key] = []
                for c1 in range(len(compute_loops_nonstiff)):
                    snap_compute_loops_nonstiff[branch_stiffness_key].append(compute_loops_nonstiff[c1])

                
                snap_loop_map_collection_nonstiff[branch_stiffness_key] = []
                for c1 in range(len(loop_map_collection_nonstiff)):
                    row_vector = []
                    for c2 in range(len(loop_map_collection_nonstiff[c1])):
                        row_vector.append(loop_map_collection_nonstiff[c1][c2])
                    snap_loop_map_collection_nonstiff[branch_stiffness_key].append(row_vector)


                snap_single_collection_stiff[branch_stiffness_key] = []
                for c1 in range(len(single_stiff_collection)):
                    row_vector = []
                    for c2 in range(len(single_stiff_collection[c1])):
                        row_vector.append(single_stiff_collection[c1][c2])
                    snap_single_collection_stiff[branch_stiffness_key].append(row_vector)
                
                
                snap_compute_loops_stiff[branch_stiffness_key] = []
                for c1 in range(len(compute_loops_stiff)):
                    snap_compute_loops_stiff[branch_stiffness_key].append(compute_loops_stiff[c1])

                
                snap_loop_map_collection_stiff[branch_stiffness_key] = []
                for c1 in range(len(loop_map_collection_stiff)):
                    row_vector = []
                    for c2 in range(len(loop_map_collection_stiff[c1])):
                        row_vector.append(loop_map_collection_stiff[c1][c2])
                    snap_loop_map_collection_stiff[branch_stiffness_key].append(row_vector)
                

                snap_nonstiff_loops[branch_stiffness_key] = []
                for c1 in range(len(nonstiff_loops)):
                    snap_nonstiff_loops[branch_stiffness_key].append(nonstiff_loops[c1])
                
            else:

                # If the layout of the branches has been found before, the
                # system matrices can be extracted from the dictionaries
                
                system_loop_map = []
                for c1 in range(len(snap_system_loopmap[branch_stiffness_key])):
                    row_vector = []
                    for c2 in range(len(snap_system_loopmap[branch_stiffness_key][c1])):
                        row_vector.append(snap_system_loopmap[branch_stiffness_key][c1][c2])
                    system_loop_map.append(row_vector)
                

                single_nonstiff_collection = []
                for c1 in range(len(snap_single_collection_nonstiff[branch_stiffness_key])):
                    row_vector = []
                    for c2 in range(len(snap_single_collection_nonstiff[branch_stiffness_key][c1])):
                        row_vector.append(snap_single_collection_nonstiff[branch_stiffness_key][c1][c2])
                    single_nonstiff_collection.append(row_vector)
                

                compute_loops_nonstiff = []
                for c1 in range(len(snap_compute_loops_nonstiff[branch_stiffness_key])):
                    compute_loops_nonstiff.append(snap_compute_loops_nonstiff[branch_stiffness_key][c1])

                
                loop_map_collection_nonstiff = []
                for c1 in range(len(snap_loop_map_collection_nonstiff[branch_stiffness_key])):
                    row_vector = []
                    for c2 in range(len(snap_loop_map_collection_nonstiff[branch_stiffness_key][c1])):
                        row_vector.append(snap_loop_map_collection_nonstiff[branch_stiffness_key][c1][c2])
                    loop_map_collection_nonstiff.append(row_vector)


                single_stiff_collection = []
                for c1 in range(len(snap_single_collection_stiff[branch_stiffness_key])):
                    row_vector = []
                    for c2 in range(len(snap_single_collection_stiff[branch_stiffness_key][c1])):
                        row_vector.append(snap_single_collection_stiff[branch_stiffness_key][c1][c2])
                    single_stiff_collection.append(row_vector)
                

                compute_loops_stiff = []
                for c1 in range(len(snap_compute_loops_stiff[branch_stiffness_key])):
                    compute_loops_stiff.append(snap_compute_loops_stiff[branch_stiffness_key][c1])

                
                loop_map_collection_stiff = []
                for c1 in range(len(snap_loop_map_collection_stiff[branch_stiffness_key])):
                    row_vector = []
                    for c2 in range(len(snap_loop_map_collection_stiff[branch_stiffness_key][c1])):
                        row_vector.append(snap_loop_map_collection_stiff[branch_stiffness_key][c1][c2])
                    loop_map_collection_stiff.append(row_vector)
                

                nonstiff_loops = []
                for c1 in range(len(snap_nonstiff_loops[branch_stiffness_key])):
                    nonstiff_loops.append(snap_nonstiff_loops[branch_stiffness_key][c1])


            # This second run of the function is to determine the new
            # branch currents because of any possible change in the state
            # of nonlinear devices.

            # Initializing the branch currents for the next set of nodal analysis
            for c1 in range(len(branch_params)):
                branch_currents[c1] = 0.0


            for c1 in range(len(branch_params)):
                kcl_branch_lookup[c1][0] = branch_events[c1]
#                if branch_events[c1]=="yes" or branch_events[c1]=="hard":
#                if branch_events[c1]=="yes" or branch_events[c1]=="hard" or branch_params[c1][-1][0][1]:
                if stiff_ratio[c1]=="yes":
                    kcl_branch_lookup[c1][3] = 1
                else:
                    if (branch_params[c1][-1][0][1]):
                        if  abs(branch_params[c1][-1][0][1]/branch_params[c1][-1][0][0]) < 0.1*dt:
                            kcl_branch_lookup[c1][3] = 1
                        else:
                            kcl_branch_lookup[c1][3] = 2
                    else:
                        kcl_branch_lookup[c1][3] = 1


                kcl_branch_lookup[c1][4] = branch_params[c1][-1][3]
                kcl_branch_lookup[c1][5] = branch_params[c1][-1][0][0]
                kcl_branch_lookup[c1][6] = branch_params[c1][-1][2]

            
            
            # Nodal analysis            
            Slv.current_continuity(kcl_branch_lookup, admittance_matrix, source_vector, abridged_node_voltage, node_voltage, \
                        kcl_node_list, node_list, shortnode_list, branch_params, nonlinear_freewheel_branches, \
                        branch_currents, "calc_currents",  t,  1)
    
            # Set the branch parameter currents equal to
            # the branch currents from the above nodal analysis.
            for c1 in range(len(branch_params)):
                branch_params[c1][-1][2] = branch_currents[c1]
            
            
            # This next part arranges the non stiff loops according
            # to their L/R ratios and attempts at isolating loops
            # with very low L/R ratio that may need to be approximated
            # as static loops. These loops will be different from
            # stiff loops as their current could ne non-negligible.
#            Slv.approximate_nonstiff_loops(branch_params, stiff_ratio, system_loop_map, branches_in_kcl_nodes, kcl_branch_map)
            single_nonstiff_collection, compute_loops_nonstiff, loop_map_collection_nonstiff, \
            single_stiff_collection, compute_loops_stiff, loop_map_collection_stiff, nonstiff_loops = \
                    Slv.compute_loop_currents_generate(branch_params, stiff_ratio, system_loop_map, branch_events, branches_in_kcl_nodes)

            
            # Compute loop currents from the branch currents
            Slv.compute_loop_currents_calc(single_nonstiff_collection, compute_loops_nonstiff, loop_map_collection_nonstiff, \
                            single_stiff_collection, compute_loops_stiff, loop_map_collection_stiff, branch_params, \
                            [curr_state_vec, next_state_vec])
            
            
            # Re-initialize the matrices A, B, and E
            # Only the nonstiff loops are solved. The stiff branches
            # are claculated by KCL.
            sys_mat_a.zeros(len(nonstiff_loops),len(nonstiff_loops))
            sys_mat_e.zeros(len(nonstiff_loops),len(nonstiff_loops))
            sys_mat_b.zeros(len(nonstiff_loops),sys_mat_b.columns)
        
            # Recalculate the system matrices for the new loops.
            Slv.recalculate_sys_matrices(system_loop_map, nonstiff_loops, branch_params, sys_mat_a, sys_mat_e, sys_mat_b, dt)

            
            # Upper triangularization of matrix E.
            Slv.mat_ode_reduce(sys_mat_e, sys_mat_a, sys_mat_b)


            # Mark the stiff loops
            stiff_loops = []
            for c1 in range(len(system_loop_map)):
                if not c1 in nonstiff_loops:
                    stiff_loops.append(c1)
                    if abs(curr_state_vec.data[c1][0])>1.0e-4:
                        curr_state_vec.data[c1][0] = 0.0
                        next_state_vec.data[c1][0] = 0.0


            # Initialize the matrices for the stiff equations
            stiff_sys_mat_a1 = matrix.Matrix(len(stiff_loops), len(stiff_loops))
            stiff_sys_mat_a2 = matrix.Matrix(len(stiff_loops),len(nonstiff_loops))
            stiff_sys_mat_e = matrix.Matrix(len(stiff_loops),len(nonstiff_loops))
            stiff_sys_mat_b = matrix.Matrix(len(stiff_loops),sys_mat_b.columns)
            
            # Set the values for the matrices for the stiff equations.
            Slv.determining_matrices_for_stiff_loops(system_loop_map, branch_params, stiff_loops, nonstiff_loops, stiff_sys_mat_a1, \
                                        stiff_sys_mat_a2, stiff_sys_mat_e, stiff_sys_mat_b)

            # Find out which are the branches and loops that have turned stiff in
            # this iteration. These branches and loops will have currents that
            # are not nelgigible but the resistances will be large.
            branches_turned_stiff = []
            loops_turned_stiff = []
            for c1 in stiff_loops:
                for c2 in range(len(system_loop_map[c1])):
                    if not system_loop_map[c1][c2]=="no":
                        if branch_events[c2]=="yes" or branch_events[c2]=="hard":
                            if stiff_ratio[c2]=="yes":
                                if c2 not in branches_turned_stiff:
                                    branches_turned_stiff.append(c2)
                                if c1 not in loops_turned_stiff:
                                    loops_turned_stiff.append(c1)
            
            
            # Reset the loops that have turned stiff in this iteration.
            for c1 in loops_turned_stiff:
                curr_state_vec.data[stiff_loops.index(c1)][0] = 0.0
                next_state_vec.data[stiff_loops.index(c1)][0] = 0.0

            
            # Determine the branches that are only in stiff loops.
#            Slv.determining_stiff_branches(system_loop_map, branches_in_stiff_loops, inductor_list, inductor_stiffness)

            # Remove the branch events
            for c1 in range(len(branch_events)):
                branch_events[c1] = "no"
    
            # Since the system has changed, generate a time_event
            # for the ODE solver to execute.
            time_event = "yes"
            
#            print "*"*100
#            print t
#            Slv.debug_loops(system_loop_map, branches_in_kcl_nodes, branch_params, branch_tags_in_loops, range(len(system_loop_map)))
#            print
#            print
            
            # End of the event driven tasks.
        
        
        # If there has been a system change or time_event
        # or if time is greater than the simulation time step.
        if ((t>=t_ode) or (time_event=="yes")):
            
            # Update the input values in u matrix
            for c1 in range(len(source_list)):
                component_objects[source_list[c1]].generate_val(source_list, sys_mat_e, sys_mat_a, sys_mat_b, sys_mat_u, t, t-t_ode_prev)
            
            # Since only the nonstiff loops are involved in the ODE
            # use a reduced order state vector and map the currents
            # to the nonstiff loops in the complete state vector.
            reduced_curr_state.zeros(len(nonstiff_loops))
            reduced_next_state.zeros(len(nonstiff_loops))
            
            for c1 in range(len(nonstiff_loops)):
                reduced_curr_state.data[c1][0] = curr_state_vec.data[nonstiff_loops[c1]][0]
                reduced_next_state.data[c1][0] = next_state_vec.data[nonstiff_loops[c1]][0]
            
            # The ODE solver.
            # Will return the x(k+1) value called
            # as next_state_vec from x(k) value
            # called curr_state_vec
            
            Slv.mat_ode(sys_mat_e, sys_mat_a, sys_mat_b, [reduced_curr_state, reduced_next_state], sys_mat_u, t-t_ode_prev, ode_var, node_list)
            
            
            # Return the currents in the reduced order state vector to the
            # nonstiff loops in the full state vector
            for c1 in range(len(nonstiff_loops)):
                curr_state_vec.data[nonstiff_loops[c1]][0] = reduced_curr_state.data[c1][0]
                next_state_vec.data[nonstiff_loops[c1]][0] = reduced_next_state.data[c1][0]


            # Calculate the currents in the stiff loops.
            for c1 in range(len(stiff_loops)-1, -1, -1):
                current_in_stiff_loop = 0
                
                for c2 in range(stiff_sys_mat_b.columns):
                    current_in_stiff_loop += stiff_sys_mat_b.data[c1][c2]*sys_mat_u.data[c2][0]
                
                for c2 in range(len(nonstiff_loops)):
                    current_in_stiff_loop -= stiff_sys_mat_a2.data[c1][c2]*curr_state_vec.data[nonstiff_loops[c2]][0]
                
                for c2 in range(len(nonstiff_loops)):
                    current_in_stiff_loop -= stiff_sys_mat_e.data[c1][c2]*ode_var[4].data[c2][0]/dt
                
                for c2 in range(c1+1, len(stiff_loops)):
                    current_in_stiff_loop -= stiff_sys_mat_a1.data[c1][c2]*curr_state_vec.data[c2][0]
                
                curr_state_vec.data[c1][0] = current_in_stiff_loop/stiff_sys_mat_a1.data[c1][c1]
                if (curr_state_vec.data[c1][0] > 0.0) and (curr_state_vec.data[c1][0] > 1.0e-4):
                    curr_state_vec.data[c1][0] = 0.0
                if (curr_state_vec.data[c1][0] < 0.0) and (curr_state_vec.data[c1][0] < -1.0e-4):
                    curr_state_vec.data[c1][0] = 0.0
                next_state_vec.data[c1][0] = curr_state_vec.data[c1][0]
            
            
            # Recalculate the branch currents from the loop currents
            for c1 in range(len(branch_params)):
                branch_params[c1][-1][4] = branch_params[c1][-1][2]
                branch_params[c1][-1][2] = 0.0
                for c2 in range(len(system_loop_map)):
                    if system_loop_map[c2][c1]=="forward" or system_loop_map[c2][c1]=="stiff_forward":
                        branch_params[c1][-1][2] += next_state_vec.data[c2][0]
    
                    elif system_loop_map[c2][c1]=="reverse" or system_loop_map[c2][c1]=="stiff_reverse":
                        branch_params[c1][-1][2] -= next_state_vec.data[c2][0]

                branch_params[c1][-1][3] = 0.0
                for c2 in range(len(branch_params[c1][-1][1])):
                    branch_params[c1][-1][3] += branch_params[c1][-1][1][c2]*sys_mat_u.data[c2][0]
            

            # This last run of KCL is to determine the currents in
            # the stiff branches.

            last_kcl_branches = []
            
            # Initializing the branch currents for the next set of nodal analysis
            for c1 in range(len(branch_params)):
                branch_currents[c1] = 0.0

            for c1 in range(len(branch_params)):
                kcl_branch_lookup[c1][0] = branch_events[c1]
#                if branch_events[c1]=="yes" or branch_events[c1]=="hard" or branch_params[c1][-1][0][1]:
                if stiff_ratio[c1]=="yes":
                    kcl_branch_lookup[c1][3] = 1
                    last_kcl_branches.append(c1)
                else:
                    if (branch_params[c1][-1][0][1]):
                        if  abs(branch_params[c1][-1][0][1]/branch_params[c1][-1][0][0]) < 0.1*dt:
                            kcl_branch_lookup[c1][3] = 1
                        else:
                            kcl_branch_lookup[c1][3] = 2
                    else:
                        kcl_branch_lookup[c1][3] = 1

                non_stiff_found = "no"
                for c2 in nonstiff_loops:
                    if not system_loop_map[c2][c1]=="no":
                        non_stiff_found = "yes"

                if non_stiff_found=="no":
                    kcl_branch_lookup[c1][3] = 1
                    if c1 not in last_kcl_branches:
                        last_kcl_branches.append(c1)

                kcl_branch_lookup[c1][4] = branch_params[c1][-1][3]
                kcl_branch_lookup[c1][5] = branch_params[c1][-1][0][0]
                kcl_branch_lookup[c1][6] = branch_params[c1][-1][2]

            
            # Nodal analysis
            Slv.current_continuity(kcl_branch_lookup, admittance_matrix, source_vector, abridged_node_voltage, node_voltage, \
                                kcl_node_list, node_list, shortnode_list, branch_params, nonlinear_freewheel_branches, branch_currents, \
                                "calc_currents",  t,  0)

            # Set the currents of the stiff branches as outputs of the KCL.
            for c1 in range(len(branch_params)):
                if stiff_ratio[c1]=="yes" or c1 in last_kcl_branches:
                    branch_params[c1][-1][2] = branch_currents[c1]
                    if branch_params[c1][-1][5]>0:
                        branch_params[c1][-1][5] -= 1
                    
                else:
                    if branch_params[c1][-1][2]>0 and branch_params[c1][-1][4]<0:
                        branch_params[c1][-1][5] += 1
                    
                    elif branch_params[c1][-1][2]<0 and branch_params[c1][-1][4]>0:
                        branch_params[c1][-1][5] += 1
                        
                    else:
                        if branch_params[c1][-1][5]>0:
                            branch_params[c1][-1][5] -= 1            
            
            
            # Update the values of objects based on x(k+1)
            for comps in component_objects.keys():
                component_objects[comps].update_val(system_loop_map, stiff_ratio, sys_mat_e, sys_mat_a, sys_mat_b, next_state_vec, \
                                    sys_mat_u, branch_params, branch_events)

            # x(k)=x(k+1) for next iteration.
            for c1 in range(system_size):
                curr_state_vec.data[c1][0] = next_state_vec.data[c1][0]
            
    
            # Save the previous time instant of ODE solver
            t_ode_prev = t
            # If the above code is running because of a
            # change in the system and not because of time
            # greater than t_ode, do not increment the t_ode.
            if time_event=="yes" and t<t_ode:
                time_event = "no"
            else:
                t_ode = t_ode + dt
            
            if (t>t_ode):
                t_ode = t_ode + dt

    
        # Store the measured values.
        if (t>=t_store):
            
            c1 = 1
            while c1<len(t_split_file):
                if t>t_split_file[c1-1] and t<t_split_file[c1]:
                    t_index = c1
                    c1 = len(t_split_file)
                else:
                    c1 += 1
            
            outputfile[t_index-1].write("%s " %str(t),)
            for c1 in range(len(meter_list)):
                outputfile[t_index-1].write("%s " %component_objects[meter_list[c1]].op_value,)
    

            if plotted_variable_list:
                for c1 in plotted_variable_list:
                    if control_file_variablestorage[c1][1].lower()=="yes" or \
                            control_file_variablestorage[c1][1].lower()=="y":
                        outputfile[t_index-1].write("%s " %control_file_variablestorage[c1][0],)
    
            outputfile[t_index-1].write("\n")
    
            t_store = t_store + dt_store

    
        # This time vector will contain all the future time events
        # generated by the control functions and the main ODE solver
        time_vector = [t_ode]
    
        # Call the control codes only if controlled elements exist.
        #if controlled_elements:
        if control_files:
            # Call the control functions in the main control programs.
            # Use the eval function to call the functions as string arguments
            for c1 in range(len(control_files)):
                eval("__control.%s(control_file_inputs, control_file_outputs, control_file_staticvars, control_file_timeevents, \
                            control_file_variablestorage, control_file_events, component_objects, c1, t, time_vector)" %control_functions[c1])

            if 1 in control_file_events:
                time_event = "yes"


        # The soonest event will be the next time instant.
        time_vector.sort()
        
        t = time_vector[0]
        if (t-t_ode_prev)>dt:
            t = t_ode_prev + dt
        
        
        # The next block was to ensure that the next time step is not too
        # close. Save the previous time instant of ODE solver
        if (t-t_ode_prev)<dt/(1/dt):
            t = t_ode_prev + t_diff
        else:
            if (t-t_ode_prev<t_diff):
                t_diff = t - t_ode_prev

    for c1 in range(len(outputfile)):
        outputfile[c1].close()

    return


simulation_iteration_collection = {}




# Create your views here.


def extract_simulation_case(request):
    """
    This function extracts the form data from a SimulationCase form
    and saves it in a SimulationCase model object.
    """
    sim_para_received = SimulationCaseForm(request.POST)
    # sim_id is always returned back and forth between view functions and
    # template files to keep track of the simulation being run. By default,
    # sim_id is -1 when a new simulation is created.
    if "sim_id" in request.POST:
        sim_id = int(request.POST["sim_id"])
        if sim_id>0:
            sim_para_model = SimulationCase.objects.get(id=sim_id)
        else:
            sim_para_model = SimulationCase()

    if "sim_state" in request.POST:
        sim_state = int(request.POST["sim_state"])

    if sim_para_received.is_valid():
        sim_parameters = sim_para_received.cleaned_data
        sim_para_model.sim_title = sim_parameters["sim_title"]
        sim_para_model.sim_descrip = sim_parameters["sim_descrip"]
        sim_para_model.sim_time_limit = sim_parameters["sim_time_limit"]
        sim_para_model.sim_time_data = sim_parameters["sim_time_data"]
        sim_para_model.sim_output_file = sim_parameters["sim_output_file"]
        sim_para_model.sim_output_slice = sim_parameters["sim_output_slice"]
        sim_para_model.sim_div_number = sim_parameters["sim_div_number"]
        sim_para_model.sim_working_directory = sim_parameters["sim_working_directory"]
        sim_para_model.save()
        sim_id = sim_para_model.id

        # simulation_form always has two parts - the simulation model and the form.
        # Only one will exist and the template file uses this property to display
        # either the form or just the parameters.
        # When parameters have been successfully saved, only the model is returned.
        simulation_form = []
        simulation_form.append(sim_para_model)
        simulation_form.append([])

    else:
        # When parameters have not been returned successfully, the form is returned
        # with errors that have been generated with custom clean methods in the model.
        simulation_form = []
        simulation_form.append([])
        simulation_form.append(sim_para_received)

    return [sim_id, sim_state, simulation_form]


def save_simulation_parameters(request):
    """
    This function saves the parameters of the simulation. If circuit
    files already have been added, it will list them or else create
    an empty circuit spreadsheet form.
    """
    # Save simulation parameters.
    sim_id, sim_state, simulation_form = extract_simulation_case(request)
    ckt_schematic_form = []
    if sim_id>0:
        sim_para_model = SimulationCase.objects.get(id=sim_id)
        # If there is no error in the form, only the model will be returned.
        if not simulation_form[1]:
            ckt_file_list = sim_para_model.circuitschematics_set.all()
            # If there are circuit files with the simulation, list them
            # and let the user add more circuits.
            if ckt_file_list:
                for ckt_file_item in ckt_file_list:
                    ckt_schematic_form.append([ckt_file_item, \
                                            CircuitSchematicsForm(instance=ckt_file_item)])
            # If no circuit files, create a blank circuit add form
            else:
                ckt_schematic_form.append([[], CircuitSchematicsForm()])
    
    return [sim_id, sim_state, simulation_form, ckt_schematic_form]


def check_control_files(sim_para_model):
    control_file_list = sim_para_model.controlfile_set.all()
    control_file_form = []
    control_errors = -1
    # List the existing circuits.
    for control_file_item in control_file_list:
        control_item_dict = model_to_dict(control_file_item)
        control_item_form = models.ControlFileForm(control_item_dict, \
                                            instance=control_file_item)
        control_full_path = os.path.join(os.sep, \
                                    sim_para_model.sim_working_directory, \
                                    control_file_item.control_file_name)

        # Try to read the file.
        try:
            check_control_file = open(control_full_path, "r")
        # If can't be read, it means file doesn't exist in the working directory.
        except:

            if control_item_form.is_valid():
                control_item_form.add_error('control_file_descrip', \
                                    'Control file could not be read. \
                                    Make sure it is in same directory as working directory above')
                control_errors = 1

            control_file_form.append([control_file_item, control_item_form])

        # If it can be read, it could be a genuine file in which case save it.
        # Or else, it may not be a .csv file, in which case raise an error.
        else:
            if len(control_file_item.control_file_name.split("."))>1 \
                    and control_file_item.control_file_name.split(".")[-1]=="py":
                control_file_form.append([control_file_item, \
                                        models.ControlFileForm(instance=control_file_item)])
            else:
                control_item_form.add_error('control_file_descrip', \
                                        'Control file must be a .py file.')
                control_file_form.append([[], control_item_form])
                control_errors = 1

    return [control_file_form, control_errors]


def save_control_func(request):
    """
    This function saves the circuit spreadsheet details and also
    runs the check whether the file exists in the working directory
    and if it is a .csv file.
    """
    if "sim_id" in request.POST:
        sim_id = int(request.POST["sim_id"])
        if sim_id>0:
            sim_para_model = SimulationCase.objects.get(id=sim_id)

    if "sim_state" in request.POST:
        sim_state = int(request.POST["sim_state"])
    
    control_file_form, control_errors = check_control_files(sim_para_model)

    # ckt_file_path in request.POST contains the circuit file name because no
    # upload takes place and only file name is obtained.
    if u'control_file' in request.POST:
        control_file = models.ControlFile()
        # Add the circuit file to the working directory path
        received_control_file_name = request.POST.getlist(u'control_file')[0]
        control_full_path = os.path.join(os.sep, \
                                sim_para_model.sim_working_directory, \
                                received_control_file_name)
        # Try to read the file.
        try:
            check_control_file = open(control_full_path, "r")
        # If can't be read, it means file doesn't exist in the working directory.
        except:
            control_form = models.ControlFileForm(request.POST)
            if control_form.is_valid():
                control_form.add_error('control_file_descrip', \
                                    'Control file could not be read. \
                                    Make sure it is in same directory as working directory above')
            control_file_form.append([[], control_form])
            control_errors = 1

        # If it can be read, it could be a genuine file in which case save it.
        # Or else, it may not be a .csv file, in which case raise an error.
        else:
            if len(received_control_file_name.split("."))>1 and \
                    received_control_file_name.split(".")[-1]=="py":
                repeated_circuit = False
                for other_control_files in sim_para_model.controlfile_set.all():
                    if received_control_file_name==other_control_files.control_file_name:
                        repeated_circuit = True
                if repeated_circuit:
                    control_form = models.ControlFileForm(request.POST)
                    control_form.add_error('control_file_descrip', \
                                    'Control file has already been added.')
                    control_file_form.append([[], control_form])
                    control_errors = 1
                else:
                    control_file.control_file_name = received_control_file_name
                    if u'control_file_descrip' in request.POST:
                        control_file.control_file_descrip = \
                                request.POST.getlist(u'control_file_descrip')[0]
                    control_file.sim_case = sim_para_model
                    control_file.save()
                    sim_para_model.save()
                    control_file_form.append([control_file, models.ControlFileForm(instance=control_file)])
            else:
                control_form = models.ControlFileForm(request.POST)
                control_form.add_error('control_file_descrip', \
                                'Control file must be a .py file.')
                control_file_form.append([[], control_form])
                control_errors = 1

    return [sim_id, sim_state, control_file_form, control_errors]


def add_control_func(request):
    """
    This function adds a blank circuit schematic form.
    """
    if "sim_id" in request.POST:
        sim_id = int(request.POST["sim_id"])
        if sim_id>0:
            sim_para_model = SimulationCase.objects.get(id=sim_id)

    if "sim_state" in request.POST:
        sim_state = int(request.POST["sim_state"])

    control_file_form, control_errors = check_control_files(sim_para_model)

    # Add a blank circuit form.
    control_file_form.append([[], models.ControlFileForm()])
    
    return [sim_id, sim_state, control_file_form, control_errors]


def check_circuit_errors(sim_para_model):
    ckt_file_list = sim_para_model.circuitschematics_set.all()
    ckt_schematic_form = []
    ckt_errors = -1
    # List the existing circuits.
    for ckt_file_item in ckt_file_list:
        ckt_item_dict = model_to_dict(ckt_file_item)
        ckt_item_form = CircuitSchematicsForm(ckt_item_dict, \
                                            instance=ckt_file_item)
        ckt_full_path = os.path.join(os.sep, \
                                    sim_para_model.sim_working_directory, \
                                    ckt_file_item.ckt_file_name)

        # Try to read the file.
        try:
            check_ckt_file = open(ckt_full_path, "r")
        # If can't be read, it means file doesn't exist in the working directory.
        except:

            if ckt_item_form.is_valid():
                ckt_item_form.add_error('ckt_file_descrip', \
                                    'Circuit spreadsheet could not be read. \
                                    Make sure it is in same directory as working directory above')
                ckt_errors = 1

            ckt_schematic_form.append([ckt_file_item, ckt_item_form])

        # If it can be read, it could be a genuine file in which case save it.
        # Or else, it may not be a .csv file, in which case raise an error.
        else:
            if len(ckt_file_item.ckt_file_name.split("."))>1 \
                    and ckt_file_item.ckt_file_name.split(".")[-1]=="csv":
                ckt_schematic_form.append([ckt_file_item, \
                                        CircuitSchematicsForm(instance=ckt_file_item)])
            else:
                ckt_item_form.add_error('ckt_file_descrip', \
                                        'Circuit schematic must be a .csv file.')
                ckt_schematic_form.append([[], ckt_item_form])
                ckt_errors = 1

    return [ckt_schematic_form, ckt_errors]


def save_circuit_schematic(request):
    """
    This function saves the circuit spreadsheet details and also
    runs the check whether the file exists in the working directory
    and if it is a .csv file.
    """
    if "sim_id" in request.POST:
        sim_id = int(request.POST["sim_id"])
        if sim_id>0:
            sim_para_model = SimulationCase.objects.get(id=sim_id)

    if "sim_state" in request.POST:
        sim_state = int(request.POST["sim_state"])

    # If the circuit schematic function is being processed,
    # only need to display simulation parameters.
    simulation_form = []
    simulation_form.append(sim_para_model)
    simulation_form.append([])
    
    ckt_schematic_form, ckt_errors = check_circuit_errors(sim_para_model)

    # ckt_file_path in request.POST contains the circuit file name because no
    # upload takes place and only file name is obtained.
    if u'ckt_file_path' in request.POST:
        ckt_file = CircuitSchematics()
        # Add the circuit file to the working directory path
        received_ckt_file_name = request.POST.getlist(u'ckt_file_path')[0]
        ckt_full_path = os.path.join(os.sep, \
                                sim_para_model.sim_working_directory, \
                                received_ckt_file_name)
        # Try to read the file.
        try:
            check_ckt_file = open(ckt_full_path, "r")
        # If can't be read, it means file doesn't exist in the working directory.
        except:
            ckt_form = CircuitSchematicsForm(request.POST)
            if ckt_form.is_valid():
                ckt_form.add_error('ckt_file_descrip', \
                                    'Circuit spreadsheet could not be read. \
                                    Make sure it is in same directory as working directory above')
            ckt_schematic_form.append([[], ckt_form])
            ckt_errors = 1

        # If it can be read, it could be a genuine file in which case save it.
        # Or else, it may not be a .csv file, in which case raise an error.
        else:
            if len(received_ckt_file_name.split("."))>1 and \
                    received_ckt_file_name.split(".")[-1]=="csv":
                repeated_circuit = False
                for other_ckt_files in sim_para_model.circuitschematics_set.all():
                    if received_ckt_file_name==other_ckt_files.ckt_file_name:
                        repeated_circuit = True
                if repeated_circuit:
                    ckt_form = CircuitSchematicsForm(request.POST)
                    ckt_form.add_error('ckt_file_descrip', \
                                    'Circuit schematic has already been added.')
                    ckt_schematic_form.append([[], ckt_form])
                    ckt_errors = 1
                else:
                    ckt_file.ckt_file_name = received_ckt_file_name
                    if u'ckt_file_descrip' in request.POST:
                        ckt_file.ckt_file_descrip = request.POST.getlist(u'ckt_file_descrip')[0]
                    ckt_file.ckt_sim_case = sim_para_model
                    ckt_file.save()
                    sim_para_model.save()
                    ckt_schematic_form.append([ckt_file, CircuitSchematicsForm(instance=ckt_file)])
            else:
                ckt_form = CircuitSchematicsForm(request.POST)
                ckt_form.add_error('ckt_file_descrip', \
                                'Circuit schematic must be a .csv file.')
                ckt_schematic_form.append([[], ckt_form])
                ckt_errors = 1

    return [sim_id, sim_state, ckt_schematic_form, ckt_errors]


def add_circuit_schematic(request):
    """
    This function adds a blank circuit schematic form.
    """
    if "sim_id" in request.POST:
        sim_id = int(request.POST["sim_id"])
        if sim_id>0:
            sim_para_model = SimulationCase.objects.get(id=sim_id)

    if "sim_state" in request.POST:
        sim_state = int(request.POST["sim_state"])

    simulation_form = []
    simulation_form.append(sim_para_model)
    simulation_form.append([])

    ckt_schematic_form, ckt_errors = check_circuit_errors(sim_para_model)

    # Add a blank circuit form.
    ckt_schematic_form.append([[], CircuitSchematicsForm()])
    
    return [sim_id, sim_state, ckt_schematic_form, ckt_errors]


def edit_simulation_parameters(request):
    """
    This function returns the user back to editing the simulation
    parameters.
    """
    if "sim_id" in request.POST:
        sim_id = int(request.POST["sim_id"])
        if sim_id>0:
            sim_para_model = SimulationCase.objects.get(id=sim_id)
        else:
            sim_para_model = SimulationCase()

    if "sim_state" in request.POST:
        sim_state = int(request.POST["sim_state"])

    simulation_form = []
    simulation_form.append([])
    simulation_form.append(SimulationCaseForm(instance=sim_para_model))

    return [sim_id, sim_state, simulation_form]


def load_simulation_parameters(request, sim_id):
    """
    This function returns the user back to editing the simulation
    parameters.
    """
    sim_para_model = SimulationCase.objects.get(id=sim_id)

    if "sim_state" in request.POST:
        sim_state = int(request.POST["sim_state"])

    simulation_form = []
    simulation_form.append([])
    simulation_form.append(SimulationCaseForm(instance=sim_para_model))

    return [sim_state, simulation_form]


def generate_component_data(sim_para_model):
    ckt_file_list = sim_para_model.circuitschematics_set.all()
    nw_input = []
    conn_ckt_mat = []
    if ckt_file_list:
        for ckt_file_item in ckt_file_list:
            nw_input.append(ckt_file_item.ckt_file_name.split(".csv")[0])
            full_file_path = os.path.join(os.sep, \
                                        sim_para_model.sim_working_directory, \
                                        ckt_file_item.ckt_file_name)
            ckt_file_object = open(full_file_path, "r")
            # Read the circuit into conn_ckt_mat
            # Also performs a scrubbing of circuit spreadsheet
            conn_ckt_mat.append(NwRdr.csv_reader(ckt_file_object))

    # Contains all the components returned by the functions below
    # needed for circuit analysis.
    circuit_analysis_components = []

    # Making a list of the type of components in the 
    # circuit.
    components_found, component_objects, ckt_error_list = \
            NwRdr.determine_circuit_components(conn_ckt_mat, nw_input)

    if not ckt_error_list:
        # Make lists of nodes and branches in the circuit.
        node_list, branch_map, node_branch_errors = \
                    NwRdr.determine_nodes_branches(conn_ckt_mat, nw_input)

        if node_branch_errors:
            ckt_error_list.extend(node_branch_errors)
        else:
            circuit_analysis_components.append(node_list)
            circuit_analysis_components.append(branch_map)

    return [components_found, component_objects, \
            circuit_analysis_components, ckt_error_list]


def process_circuit_schematics(sim_para_model):
    ckt_file_list = sim_para_model.circuitschematics_set.all()
    nw_input = []
    conn_ckt_mat = []
    if ckt_file_list:
        for ckt_file_item in ckt_file_list:
            nw_input.append(ckt_file_item.ckt_file_name.split(".csv")[0])
            full_file_path = os.path.join(os.sep, \
                                        sim_para_model.sim_working_directory, \
                                        ckt_file_item.ckt_file_name)
            ckt_file_object = open(full_file_path, "r")
            # Read the circuit into conn_ckt_mat
            # Also performs a scrubbing of circuit spreadsheet
            conn_ckt_mat.append(NwRdr.csv_reader(ckt_file_object))

    # Contains all the components returned by the functions below
    # needed for circuit analysis.
    circuit_analysis_components = []

    # Making a list of the type of components in the 
    # circuit.
    components_found, component_objects, ckt_error_list = \
            NwRdr.determine_circuit_components(conn_ckt_mat, nw_input)
    
    if not ckt_error_list:
        all_components = sim_para_model.circuitcomponents_set.all()
        for comp_types in components_found.keys():
            # Take every type of component found
            # item -> resistor, inductor etc
            for c1 in range(len(components_found[comp_types])):
                # Each component type will be occurring
                # multiple times. Iterate through every find.
                # The list corresponding to each component is
                # the unique cell position in the spreadsheet
                check_comp_exists = all_components.filter(comp_type=comp_types).\
                        filter(comp_tag=components_found[comp_types][c1][1])
                if check_comp_exists and len(check_comp_exists)==1:
                    old_comp_object = check_comp_exists[0]
                    old_comp_object.comp_number = c1 + 1
                    old_comp_object.comp_pos_3D = components_found[comp_types][c1][0]
                    old_comp_object.comp_pos = NwRdr.csv_element_2D(\
                                        NwRdr.csv_tuple(components_found[comp_types][c1][0])[1:])
                    sheet_number = NwRdr.csv_tuple(components_found[comp_types][c1][0])[0]
                    old_comp_object.comp_sheet = sheet_number
                    old_comp_object.sheet_name = nw_input[sheet_number] + ".csv"
                    old_comp_object.sim_case = sim_para_model
                    old_comp_object.save()
                    sim_para_model.save()
                else:
                    new_comp_object = CircuitComponents()
                    new_comp_object.comp_type = comp_types
                    new_comp_object.comp_number = c1 + 1
                    new_comp_object.comp_pos_3D = components_found[comp_types][c1][0]
                    new_comp_object.comp_pos = NwRdr.csv_element_2D(\
                                        NwRdr.csv_tuple(components_found[comp_types][c1][0])[1:])
                    sheet_number = NwRdr.csv_tuple(components_found[comp_types][c1][0])[0]
                    new_comp_object.comp_sheet = sheet_number
                    new_comp_object.sheet_name = nw_input[sheet_number] + ".csv"
                    new_comp_object.comp_tag = components_found[comp_types][c1][1]
                    new_comp_object.sim_case = sim_para_model
                    new_comp_object.save()
                    sim_para_model.save()
        
        for c2 in range(len(all_components)-1, -1, -1):
            comp_exist = all_components[c2]
            comp_found = False
            c1 = 0
            while c1<len(components_found[comp_exist.comp_type]) and comp_found==False:
                if components_found[comp_exist.comp_type][c1][1]==comp_exist.comp_tag:
                    comp_found = True
                c1 += 1
            
            if comp_found==False:
                comp_exist.delete()
                sim_para_model.save()


    if not ckt_error_list:
        # Make lists of nodes and branches in the circuit.
        node_list, branch_map, node_branch_errors = \
                    NwRdr.determine_nodes_branches(conn_ckt_mat, nw_input)

        if node_branch_errors:
            ckt_error_list.extend(node_branch_errors)
        else:
            circuit_analysis_components.append(node_list)
            circuit_analysis_components.append(branch_map)

    return [components_found, component_objects, \
            circuit_analysis_components, ckt_error_list]


def generating_control_comp_list(sim_para_model, config_control_file, control_form_type):
    control_component_list = []
    try:
        control_input_list = \
                config_control_file.controlinputs_set.all()
    except:
        control_input_list = []
    input_component_list = []
    if control_input_list:
        for input_item in control_input_list:
            input_component_list.append([input_item, []])
    if control_form_type==1:
        input_component_list.append([[], models.ControlInputsForm()])
    control_component_list.append(input_component_list)

    try:
        control_output_list = \
                config_control_file.controloutputs_set.all()
    except:
        control_output_list = []
    output_component_list = []
    if control_output_list:
        for output_item in control_output_list:
            output_component_list.append([output_item, []])
    if control_form_type==2:
        output_component_list.append([[], models.ControlOutputsForm()])
    control_component_list.append(output_component_list)

    try:
        control_staticvar_list = \
                config_control_file.controlstaticvariable_set.all()
    except:
        control_staticvar_list = []
    staticvar_component_list = []
    if control_staticvar_list:
        for staticvar_item in control_staticvar_list:
            staticvar_component_list.append([staticvar_item, []])
    if control_form_type==3:
        staticvar_component_list.append([[], models.ControlStaticVariableForm()])
    control_component_list.append(staticvar_component_list)

    try:
        control_timeevent_list = \
                config_control_file.controltimeevent_set.all()
    except:
        control_timeevent_list = []
    timeevent_component_list = []
    if control_timeevent_list:
        for timeevent_item in control_timeevent_list:
            timeevent_component_list.append([timeevent_item, []])
    if control_form_type==4:
        timeevent_component_list.append([[], models.ControlTimeEventForm()])
    control_component_list.append(timeevent_component_list)

    try:
        control_varstore_list = \
                sim_para_model.controlvariablestorage_set.all().\
                filter(control_file_name=config_control_file.control_file_name)
    except:
        control_varstore_list = []
    varstore_component_list = []
    if control_varstore_list:
        for varstore_item in control_varstore_list:
            varstore_component_list.append([varstore_item, []])
    if control_form_type==5:
        varstore_component_list.append([[], models.ControlVariableStorageForm()])
    control_component_list.append(varstore_component_list)

    return control_component_list


def index(request):
    return render(request, "index.html")


def new_simulation(request):
    if not request.method == "POST":
        # When starting a new simulation, set the sim_id to -1.
        sim_id = -1
        sim_state = 1
        simulation_form = []
        simulation_form.append([])
        simulation_form.append(SimulationCaseForm())
        ckt_schematic_form = []
        return render(request,
                "edit_simulation.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'simulation_form' : simulation_form})

    else:
        print(request.POST)
        print
        print
        
        if "save_ckt_schematic" in request.POST and \
                request.POST["save_ckt_schematic"]=="Save circuit file":
            sim_id, sim_state, ckt_schematic_form, ckt_errors = save_circuit_schematic(request)
            ckt_error_list = []
            if sim_state==1:
                sim_state = 2
            return render(request,
                "edit_circuit_schematic.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'ckt_schematic_form' : ckt_schematic_form,
                'ckt_errors' : ckt_errors,
                'ckt_error_list' : ckt_error_list})
            

        elif "add_ckt_schematic" in request.POST and \
                request.POST["add_ckt_schematic"]=="Add circuit file":
            sim_id, sim_state, ckt_schematic_form, ckt_errors = add_circuit_schematic(request)
            ckt_error_list = []
            if sim_state==1:
                sim_state = 2
            return render(request,
                "edit_circuit_schematic.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'ckt_schematic_form' : ckt_schematic_form,
                'ckt_errors' : ckt_errors,
                'ckt_error_list' : ckt_error_list})

        elif "save_sim_param" in request.POST and \
                request.POST["save_sim_param"]=="Save Simulation Parameters":
            sim_id, sim_state, simulation_form, ckt_schematic_form = \
                    save_simulation_parameters(request)
            if sim_state==1:
                sim_state = 2
            return render(request,
                    "edit_simulation.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'simulation_form' : simulation_form})

        elif "edit_sim_param" in request.POST and \
                request.POST["edit_sim_param"]=="Edit Simulation Parameters":
            sim_id, sim_state, simulation_form = edit_simulation_parameters(request)
            return render(request,
                    "edit_simulation.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'simulation_form' : simulation_form})
        
        elif "main_page" in request.POST and request.POST["main_page"]=="Back to main page":
            error_codes = []
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    ckt_file_list = sim_para_model.circuitschematics_set.all()
                    if ckt_file_list:
                        for ckt_file_item in ckt_file_list:
                            full_file_path = os.path.join(os.sep, \
                                                    sim_para_model.sim_working_directory, \
                                                    ckt_file_item.ckt_file_name)
                            try:
                                try_read_file = open(full_file_path, "r")
                            except:
                                if 1 not in error_codes:
                                    error_codes.append(1)

                else:
                    sim_para_model = SimulationCase()
            
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])

            return render(request,
                    "new_simulation.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'error_codes' : error_codes})

        elif "process_ckt_schematics" in request.POST and \
                request.POST["process_ckt_schematics"]=="Process circuit schematics":
            error_codes = []
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id<0:
                    if 2 not in error_codes:
                        error_codes.append(2)
                    sim_para_model = SimulationCase()
                else:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)

            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])

            ckt_file_list = sim_para_model.circuitschematics_set.all()
            ckt_schematic_form = []
            if ckt_file_list:
                for ckt_file_item in ckt_file_list:
                    ckt_schematic_form.append([ckt_file_item, \
                                            CircuitSchematicsForm(instance=ckt_file_item)])

            components_found, component_objects, \
                    circuit_analysis_components, ckt_error_list = \
                    process_circuit_schematics(sim_para_model)

            try:
                meter_list = sim_para_model.metercomponents_set.all()
            except:
                meter_list = []

            for comp_items in component_objects.keys():
                if component_objects[comp_items].is_meter=="yes":
                    meter_found = False
                    if meter_list:
                        check_meter = meter_list.\
                                filter(comp_type=component_objects[comp_items].type).\
                                filter(comp_tag=component_objects[comp_items].tag)
                        if check_meter and len(check_meter)==1:
                            old_meter_item = check_meter[0]
                            old_meter_item.ckt_file_name = \
                                    component_objects[comp_items].sheet_name
                            old_meter_item.comp_pos_3D = \
                                    component_objects[comp_items].pos_3D
                            old_meter_item.save()
                            sim_para_model.save()
                            meter_found = True
                        else:
                            meter_found = False
                    if not meter_found:
                        new_meter_item = models.MeterComponents()
                        new_meter_item.comp_type = component_objects[comp_items].type
                        new_meter_item.comp_tag = component_objects[comp_items].tag
                        new_meter_item.comp_pos_3D = component_objects[comp_items].pos_3D
                        new_meter_item.ckt_file_name = component_objects[comp_items].sheet_name
                        new_meter_item.comp_name = component_objects[comp_items].type + \
                                "_" + component_objects[comp_items].tag
                        new_meter_item.sim_case = sim_para_model
                        new_meter_item.save()
                        sim_para_model.save()

            try:
                control_comp_list = sim_para_model.controllablecomponents_set.all()
            except:
                control_comp_list = []

            for comp_items in component_objects.keys():
                if component_objects[comp_items].has_control=="yes":
                    controllable_comp_found = False
                    if control_comp_list:
                        check_control_comp = control_comp_list.\
                                filter(comp_type=component_objects[comp_items].type).\
                                filter(comp_tag=component_objects[comp_items].tag)
                        if check_control_comp and len(check_control_comp)==1:
                            old_control_item = check_control_comp[0]
                            old_control_item.ckt_file_name = \
                                    component_objects[comp_items].sheet_name
                            old_control_item.comp_pos_3D = \
                                    component_objects[comp_items].pos_3D
                            old_control_item.save()
                            sim_para_model.save()
                            controllable_comp_found = True
                        else:
                            controllable_comp_found = False
                    if not controllable_comp_found:
                        new_control_item = models.ControllableComponents()
                        new_control_item.comp_type = component_objects[comp_items].type
                        new_control_item.comp_tag = component_objects[comp_items].tag
                        new_control_item.comp_pos_3D = component_objects[comp_items].pos_3D
                        new_control_item.ckt_file_name = component_objects[comp_items].sheet_name
                        new_control_item.comp_name = component_objects[comp_items].type + \
                                "_" + component_objects[comp_items].tag
                        new_control_item.control_tag = component_objects[comp_items].control_tag
                        new_control_item.sim_case = sim_para_model
                        new_control_item.save()
                        sim_para_model.save()

            if not ckt_error_list:
                sim_state = 3
                ckt_errors = 0
            else:
                ckt_errors = 1

            return render(request,
                "edit_circuit_schematic.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'ckt_schematic_form' : ckt_schematic_form,
                'ckt_errors' : ckt_errors,
                'ckt_error_list' : ckt_error_list})

        elif "edit_ckt_parameters" in request.POST and \
                (request.POST["edit_ckt_parameters"]=="Edit circuit parameters" or \
                request.POST["edit_ckt_parameters"]=="Back to circuit list"):
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])

            sim_para_model = SimulationCase.objects.get(id=sim_id)
            components_found, component_objects, \
                    circuit_analysis_components, ckt_error_list = \
                    generate_component_data(sim_para_model)

            all_ckt_component_list = []
            if not ckt_error_list:
                ckt_file_list = sim_para_model.circuitschematics_set.all()
                branch_map = circuit_analysis_components[1]
                for comp_items in component_objects.keys():
                    component_objects[comp_items].\
                            create_form_values(sim_para_model, ckt_file_list, branch_map)

                for ckt_file_item in ckt_file_list:
                    all_ckt_component_list.append(ckt_file_item)

            sim_state = 4

            return render(request,
                "main_circuit_components.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'ckt_schematics_update' : all_ckt_component_list,
                'ckt_error_list' : ckt_error_list})

        elif "view_output" in request.POST and request.POST["view_output"]=="View output":
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            
            ckt_error_list = []
            circuit_analysis_components = []
            if sim_id>0:
                sim_para_model = SimulationCase.objects.get(id=sim_id)
                ckt_file_list = sim_para_model.circuitschematics_set.all()
                nw_input = []
                conn_ckt_mat = []
                for ckt_file_item in ckt_file_list:
                    ckt_full_path = os.path.join(os.sep, \
                                    sim_para_model.sim_working_directory, \
                                    ckt_file_item.ckt_file_name)

                    # Try to read the file.
                    try:
                        ckt_file_object = open(ckt_full_path, "r")
                    except:
                        ckt_error_list.append("Circuit schematic spreadsheet"+\
                            ckt_file_item.ckt_file_name+"cannot be read.")
                    else:
                        nw_input.append(ckt_file_item.ckt_file_name.split(".csv")[0])
                        # Read the circuit into conn_ckt_mat
                        # Also performs a scrubbing of circuit spreadsheet
                        conn_ckt_mat.append(NwRdr.csv_reader(ckt_file_object))

                if not ckt_error_list:
                    # Making a list of the type of components in the 
                    # circuit.
                    components_found, component_objects, network_error_list = \
                            NwRdr.determine_circuit_components(conn_ckt_mat, nw_input)

                    ckt_error_list.extend(network_error_list)

                if not ckt_error_list:
                    # Make lists of nodes and branches in the circuit.
                    node_list, branch_map, node_branch_errors = \
                                NwRdr.determine_nodes_branches(conn_ckt_mat, nw_input)

                    if node_branch_errors:
                        ckt_error_list.extend(node_branch_errors)
                    else:
                        circuit_analysis_components.append(node_list)
                        circuit_analysis_components.append(branch_map)

                if not ckt_error_list:
                    for comp_keys in component_objects.keys():
                        comp_error = component_objects[comp_keys].\
                            pre_run_check(ckt_file_item, circuit_analysis_components[1])
                        if comp_error:
                            ckt_error_list.extend(comp_error)

                if not ckt_error_list:
                    if "sim"+str(sim_para_model.id) in simulation_iteration_collection.keys():
                        run_state = 1
                    else:
                        run_state = 0

                plot_form_list = []
                for ckt_plots in sim_para_model.circuitplot_set.all():
                    prev_plot_list = []
                    for plot_items in ckt_plots.circuitwaveforms_set.all():
                        prev_plot_list.append([plot_items, 0])
                    plot_form_list.append([[ckt_plots, 0], prev_plot_list, 0])

                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'run_state' : run_state,
                    'plot_id' : -1,
                    'ckt_error_list' : ckt_error_list,
                    'plot_form_list' : plot_form_list})

        elif "run_simulation" in request.POST and request.POST["run_simulation"]=="Run":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)

                    ckt_file_list = sim_para_model.circuitschematics_set.all()
                    nw_input = []
                    conn_ckt_mat = []
                    ckt_error_list = []
                    for ckt_file_item in ckt_file_list:
                        ckt_full_path = os.path.join(os.sep, \
                                        sim_para_model.sim_working_directory, \
                                        ckt_file_item.ckt_file_name)

                        # Try to read the file.
                        try:
                            ckt_file_object = open(ckt_full_path, "r")
                        except:
                            ckt_error_list.append("Circuit schematic spreadsheet"+\
                                ckt_file_item.ckt_file_name+"cannot be read.")
                        else:
                            nw_input.append(ckt_file_item.ckt_file_name.split(".csv")[0])
                            # Read the circuit into conn_ckt_mat
                            # Also performs a scrubbing of circuit spreadsheet
                            conn_ckt_mat.append(NwRdr.csv_reader(ckt_file_object))

                    if not ckt_error_list:
                        # Making a list of the type of components in the 
                        # circuit.
                        components_found, component_objects, network_error_list = \
                                NwRdr.determine_circuit_components(conn_ckt_mat, nw_input)
                        sim_components = [components_found, \
                                component_objects]
                        ckt_error_list.extend(network_error_list)

                    if not ckt_error_list:
                        # Make lists of nodes and branches in the circuit.
                        node_list, branch_map, node_branch_errors = \
                                    NwRdr.determine_nodes_branches(conn_ckt_mat, nw_input)

                        if node_branch_errors:
                            ckt_error_list.extend(node_branch_errors)
                        else:
                            circuit_analysis_components = [node_list, branch_map]

                    if not ckt_error_list:
                        for comp_keys in component_objects.keys():
                            component_objects[comp_keys].\
                                assign_parameters(ckt_file_list)

                    if not ckt_error_list:
                        for comp_keys in component_objects.keys():
                            comp_error = component_objects[comp_keys].\
                                pre_run_check(ckt_file_item, circuit_analysis_components[1])
                            if comp_error:
                                ckt_error_list.extend(comp_error)

                    if not ckt_error_list:
                        synthesized_ckt_comps = prepare_simulation_objects(sim_para_model, \
                                    sim_components, \
                                    circuit_analysis_components, \
                                    conn_ckt_mat)
                        meter_list = synthesized_ckt_comps["meter_list"]
                        
                        try:
                            all_circuit_plotlines = sim_para_model.plotlines_set.all()
                        except:
                            new_plotline = models.PlotLines()
                            new_plotline.line_name = \
                                            component_objects[meter_list[0]].type + "_" + \
                                                    component_objects[meter_list[0]].tag
                            new_plotline.line_type = "M"
                            new_plotline.line_pos = 1
                            new_plotline.sim_case = sim_para_model
                            new_plotline.save()
                            sim_para_model.save()

                        for meter_item in meter_list:
                            meter_waveform = \
                                    component_objects[meter_item].type + "_" + \
                                    component_objects[meter_item].tag
                            check_plotline = \
                                    all_circuit_plotlines.filter(line_type="M").\
                                    filter(line_name=meter_waveform)

                            if check_plotline and len(check_plotline)==1:
                                old_plotline = check_plotline[0]
                                old_plotline.line_pos = meter_list.index(meter_item) + 1
                                old_plotline.save()
                                sim_para_model.save()
                            else:
                                new_plotline = models.PlotLines()
                                new_plotline.line_name = \
                                        component_objects[meter_item].type + "_" + \
                                        component_objects[meter_item].tag
                                new_plotline.line_type = "M"
                                new_plotline.line_pos = meter_list.index(meter_item) + 1
                                new_plotline.sim_case = sim_para_model
                                new_plotline.save()
                                sim_para_model.save()

                            meter_circuit_plotlines = all_circuit_plotlines.filter(line_type="M")
                            for c1 in range(len(meter_circuit_plotlines)-1, -1, -1):
                                meter_extra = True
                                for c2 in range(len(meter_list)):
                                    meter_plot_name = component_objects[meter_list[c2]].type + "_" + \
                                            component_objects[meter_list[c2]].tag
                                    if meter_circuit_plotlines[c1].line_name==meter_plot_name:
                                        meter_extra = False
                                if meter_extra:
                                    meter_circuit_plotlines[c1].delete()
                                    sim_para_model.save()

                        
                        plotted_variable_list = synthesized_ckt_comps["plotted_variable_list"]
                        for varstore_item in plotted_variable_list:
                            check_plotline = \
                                    all_circuit_plotlines.filter(line_type="V").\
                                    filter(line_name=varstore_item)

                            if check_plotline and len(check_plotline)==1:
                                old_plotline = check_plotline[0]
                                old_plotline.line_pos = len(meter_list) + \
                                        plotted_variable_list.index(varstore_item) + 1
                                old_plotline.save()
                                sim_para_model.save()
                            else:
                                new_plotline = models.PlotLines()
                                new_plotline.line_name = varstore_item
                                new_plotline.line_type = "V"
                                new_plotline.line_pos = len(meter_list) + \
                                        plotted_variable_list.index(varstore_item) + 1
                                new_plotline.sim_case = sim_para_model
                                new_plotline.save()
                                sim_para_model.save()
                        
                        for check_plots in sim_para_model.plotlines_set.all():
                            print(check_plots.line_name)
                            print(check_plots.line_pos)
                            print(check_plots.line_type)
                            print
                            print

                        components_in_branch = synthesized_ckt_comps["components_in_branch"]
                        comp_types_with_resistance = ["Resistor", \
                                    "Voltmeter", \
                                    "Diode", \
                                    "Switch", \
                                    "VariableResistor"]
                        for branch_check in components_in_branch:
                            resistor_found = False
                            for branch_comps in branch_check:
                                if component_objects[branch_comps].type in comp_types_with_resistance:
                                    resistor_found = True
                            if not resistor_found:
                                branch_error = "Branch with "
                                for branch_comps in branch_check:
                                    branch_error += component_objects[branch_comps].type + \
                                            "_" + component_objects[branch_comps].tag
                                    if len(branch_check)>1:
                                        branch_error += ", "
                                branch_error += " does not have a resistor."
                                ckt_error_list.append(branch_error)

                if not ckt_error_list:
                    if "sim"+str(sim_para_model.id) not in simulation_iteration_collection.keys():
                        simulator_loop = multiprocessing.Process(target=simulation_iterations, \
                                kwargs={'sim_id':sim_id, \
                                    'synthesized_ckt_comps':synthesized_ckt_comps, })
                        simulator_loop.start()
                        simulation_iteration_collection["sim"+str(sim_para_model.id)] = \
                                simulator_loop
                    run_state = 1
                else:
                    run_state = 0

                plot_form_list = []
                for ckt_plots in sim_para_model.circuitplot_set.all():
                    prev_plot_list = []
                    for plot_items in ckt_plots.circuitwaveforms_set.all():
                        prev_plot_list.append([plot_items, 0])
                    plot_form_list.append([[ckt_plots, 0], prev_plot_list, 0])

                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'run_state' : run_state,
                    'plot_id' : -1,
                    'ckt_error_list' : ckt_error_list,
                    'plot_form_list' : plot_form_list})

        elif "stop_simulation" in request.POST and request.POST["stop_simulation"]=="Stop":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)

                    if "sim"+str(sim_para_model.id) in simulation_iteration_collection.keys():
                        simulator_loop = \
                                simulation_iteration_collection["sim"+str(sim_para_model.id)]
                        simulator_loop.terminate()
                        del simulation_iteration_collection["sim"+str(sim_para_model.id)]
                
                ckt_error_list = []

                run_state = 2

                plot_form_list = []
                for ckt_plots in sim_para_model.circuitplot_set.all():
                    prev_plot_list = []
                    for plot_items in ckt_plots.circuitwaveforms_set.all():
                        prev_plot_list.append([plot_items, 0])
                    plot_form_list.append([[ckt_plots, 0], prev_plot_list, 0])

                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'run_state' : run_state,
                    'plot_id' : -1,
                    'ckt_error_list' : ckt_error_list,
                    'plot_form_list' : plot_form_list})

        elif "add_plot" in request.POST and request.POST["add_plot"]=="Add plot":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)

                    if "sim"+str(sim_para_model.id) in simulation_iteration_collection.keys():
                        run_state = 1
                    else:
                        run_state = 0

                    plot_form_list = []
                    for ckt_plots in sim_para_model.circuitplot_set.all():
                        prev_plot_list = []
                        for plot_items in ckt_plots.circuitwaveforms_set.all():
                            prev_plot_list.append([plot_items, 0])
                        plot_form_list.append([[ckt_plots, 0], prev_plot_list, 0])

                    new_circuit_form = models.CircuitPlotForm()
                    plot_form_list.append([[new_circuit_form, 1], [], 1])

                ckt_error_list = []

                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'run_state' : run_state,
                    'plot_id' : -1,
                    'ckt_error_list' : ckt_error_list,
                    'plot_form_list' : plot_form_list})

        elif "start_plot" in request.POST and request.POST["start_plot"]=="Start plot":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)

                    if "sim"+str(sim_para_model.id) in simulation_iteration_collection.keys():
                        run_state = 1
                    else:
                        run_state = 0

                    plot_form_list = []

                    received_circuit_form = models.CircuitPlotForm(request.POST)
                    if received_circuit_form.is_valid():
                        new_circuit_plot = models.CircuitPlot()
                        received_circuit_data = received_circuit_form.cleaned_data
                        new_circuit_plot.plot_title = received_circuit_data["plot_title"]
                        new_circuit_plot.sim_case = sim_para_model
                        new_circuit_plot.save()
                        sim_para_model.save()

                    for ckt_plots in sim_para_model.circuitplot_set.all():
                        if not ckt_plots.id==new_circuit_plot.id:
                            prev_plot_list = []
                            for plot_items in ckt_plots.circuitwaveforms_set.all():
                                prev_plot_list.append([plot_items, 0])
                            plot_form_list.append([[ckt_plots, 0], prev_plot_list, 0])
                    
                    plot_form_list.append([[new_circuit_plot, 0], [], 1])


                ckt_error_list = []

                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'run_state' : run_state,
                    'plot_id' : new_circuit_plot.id,
                    'ckt_error_list' : ckt_error_list,
                    'plot_form_list' : plot_form_list})

        elif "add_waveform" in request.POST and request.POST["add_waveform"]=="Add waveform":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)

                    if "sim"+str(sim_para_model.id) in simulation_iteration_collection.keys():
                        run_state = 1
                    else:
                        run_state = 0

                    plot_form_list = []

                    all_waveforms_sim = sim_para_model.plotlines_set.all()

                    if "plot_id" in request.POST:
                        plot_id = int(request.POST["plot_id"])
                        new_circuit_plot = models.CircuitPlot.objects.get(id=plot_id)

                    for ckt_plots in sim_para_model.circuitplot_set.all():
                        if not ckt_plots.id==new_circuit_plot.id:
                            prev_plot_list = []
                            for plot_items in ckt_plots.circuitwaveforms_set.all():
                                prev_plot_list.append([plot_items, 0])
                            plot_form_list.append([[ckt_plots, 0], prev_plot_list, 0])

                    prev_plot_list = []
                    for plot_items in new_circuit_plot.circuitwaveforms_set.all():
                        prev_plot_list.append([plot_items, 0])
                    new_waveform_form = models.CircuitWaveformsForm()
                    prev_plot_list.append([new_waveform_form, 1])
                    plot_form_list.append([[new_circuit_plot, 0], prev_plot_list, 1])

                ckt_error_list = []

                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'run_state' : run_state,
                    'plot_id' : plot_id,
                    'ckt_error_list' : ckt_error_list,
                    'plot_form_list' : plot_form_list,
                    'all_waveforms_sim' : all_waveforms_sim})

        elif "save_waveform" in request.POST and request.POST["save_waveform"]=="Save waveform":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)

                    if "sim"+str(sim_para_model.id) in simulation_iteration_collection.keys():
                        run_state = 1
                    else:
                        run_state = 0

                    plot_form_list = []

                    all_waveforms_sim = sim_para_model.plotlines_set.all()
                    if "plot_id" in request.POST:
                        plot_id = int(request.POST["plot_id"])
                        new_circuit_plot = models.CircuitPlot.objects.get(id=plot_id)
                    
                    new_waveform = models.CircuitWaveforms()
                    new_waveform.circuit_plot = new_circuit_plot
                    received_waveform_form = models.CircuitWaveformsForm(request.POST)
                    if received_waveform_form.is_valid():
                        received_waveform_data = received_waveform_form.cleaned_data
                        new_waveform.waveform_legend = received_waveform_data["waveform_legend"]
                        new_waveform.waveform_scale = int(received_waveform_data["waveform_scale"])
                        new_waveform.save()
                    if "waveform_source" in request.POST:
                        new_plotline = request.POST["waveform_source"]
                        new_plotline_object = sim_para_model.plotlines_set.\
                                filter(line_type="M").\
                                filter(line_name=new_plotline)
                        if new_plotline_object and len(new_plotline_object)==1:
                            new_waveform.waveform.add(new_plotline_object[0])
                        else:
                            other_plotline_object = sim_para_model.plotlines_set.\
                                filter(line_type="V").\
                                filter(line_name=new_plotline)
                            if other_plotline_object and len(other_plotline_object)==1:
                                new_waveform.waveform.add(other_plotline_object[0])
                    new_waveform.save()
                    new_circuit_plot.save()
                    sim_para_model.save()

                    for ckt_plots in sim_para_model.circuitplot_set.all():
                        if not ckt_plots.id==new_circuit_plot.id:
                            prev_plot_list = []
                            for plot_items in ckt_plots.circuitwaveforms_set.all():
                                prev_plot_list.append([plot_items, 0])
                            plot_form_list.append([[ckt_plots, 0], prev_plot_list, 0])

                    prev_plot_list = []
                    for plot_items in new_circuit_plot.circuitwaveforms_set.all():
                        prev_plot_list.append([plot_items, 0])
                    plot_form_list.append([[new_circuit_plot, 0], prev_plot_list, 1])

                ckt_error_list = []

                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'run_state' : run_state,
                    'plot_id' : plot_id,
                    'ckt_error_list' : ckt_error_list,
                    'plot_form_list' : plot_form_list,
                    'all_waveforms_sim' : all_waveforms_sim})

        elif "save_plot" in request.POST and request.POST["save_plot"]=="Save plot":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)

                    if "sim"+str(sim_para_model.id) in simulation_iteration_collection.keys():
                        run_state = 1
                    else:
                        run_state = 0

                    plot_form_list = []
                    for ckt_plots in sim_para_model.circuitplot_set.all():
                        prev_plot_list = []
                        for plot_items in ckt_plots.circuitwaveforms_set.all():
                            prev_plot_list.append([plot_items, 0])
                        plot_form_list.append([[ckt_plots, 0], prev_plot_list, 0])

                ckt_error_list = []

                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'run_state' : run_state,
                    'plot_id' : -1,
                    'ckt_error_list' : ckt_error_list,
                    'plot_form_list' : plot_form_list})

        if "save_control" in request.POST and \
                request.POST["save_control"]=="Save control file":
            sim_id, sim_state, control_file_form, ckt_errors = save_control_func(request)
            ckt_error_list = []
            return render(request,
                "add_control_files.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'control_file_form' : control_file_form,
                'ckt_errors' : ckt_errors,
                'ckt_error_list' : ckt_error_list})
            

        elif "add_control" in request.POST and \
                (request.POST["add_control"]=="Add control file" or \
                request.POST["add_control"]=="Edit control"):
            sim_id, sim_state, control_file_form, ckt_errors = add_control_func(request)
            ckt_error_list = []
            return render(request,
                "add_control_files.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'control_file_form' : control_file_form,
                'ckt_errors' : ckt_errors,
                'ckt_error_list' : ckt_error_list})

        elif "add_control_input" in request.POST and \
                request.POST["add_control_input"]=="Add control input":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "control_id" in request.POST:
                control_id = request.POST["control_id"]
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    config_control_file = sim_para_model.controlfile_set.get(id=control_id)

                    try:
                        meter_list = sim_para_model.metercomponents_set.all()
                    except:
                        meter_list = []

                    meter_name_list = []
                    for meter_item in meter_list:
                        meter_name_list.append(meter_item)

                    control_component_list = generating_control_comp_list(\
                                    sim_para_model, config_control_file, 1)

                    return render(request,
                        "config_control_files.html",
                        {'sim_id' : sim_id,
                        'sim_state' : sim_state,
                        'control_id' : control_id,
                        'control_component_list' : control_component_list,
                        'meter_list' : meter_name_list})

        elif "add_control_output" in request.POST and \
                request.POST["add_control_output"]=="Add control output":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "control_id" in request.POST:
                control_id = request.POST["control_id"]
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    config_control_file = sim_para_model.controlfile_set.get(id=control_id)

                    try:
                        control_model_list = sim_para_model.controllablecomponents_set.all()
                    except:
                        control_model_list = []

                    controlmodel_name_list = []
                    for control_model_item in control_model_list:
                        controlmodel_name_list.append(control_model_item)

                    control_component_list = generating_control_comp_list(\
                                    sim_para_model, config_control_file, 2)

                    return render(request,
                        "config_control_files.html",
                        {'sim_id' : sim_id,
                        'sim_state' : sim_state,
                        'control_id' : control_id,
                        'control_component_list' : control_component_list,
                        'control_model_list' : controlmodel_name_list})

        elif "add_control_static" in request.POST and \
                request.POST["add_control_static"]=="Add static variable":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "control_id" in request.POST:
                control_id = request.POST["control_id"]
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    config_control_file = sim_para_model.controlfile_set.get(id=control_id)

                    control_component_list = generating_control_comp_list(\
                                    sim_para_model, config_control_file, 3)

                    return render(request,
                        "config_control_files.html",
                        {'sim_id' : sim_id,
                        'sim_state' : sim_state,
                        'control_id' : control_id,
                        'control_component_list' : control_component_list})

        elif "add_control_timeevent" in request.POST and \
                request.POST["add_control_timeevent"]=="Add time event":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "control_id" in request.POST:
                control_id = request.POST["control_id"]
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    config_control_file = sim_para_model.controlfile_set.get(id=control_id)

                    control_component_list = generating_control_comp_list(\
                                    sim_para_model, config_control_file, 4)

                    return render(request,
                        "config_control_files.html",
                        {'sim_id' : sim_id,
                        'sim_state' : sim_state,
                        'control_id' : control_id,
                        'control_component_list' : control_component_list})

        elif "add_control_varstore" in request.POST and \
                request.POST["add_control_varstore"]=="Add variable storage":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "control_id" in request.POST:
                control_id = request.POST["control_id"]
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    config_control_file = sim_para_model.controlfile_set.get(id=control_id)

                    control_component_list = generating_control_comp_list(\
                                    sim_para_model, config_control_file, 5)

                    return render(request,
                        "config_control_files.html",
                        {'sim_id' : sim_id,
                        'sim_state' : sim_state,
                        'control_id' : control_id,
                        'control_component_list' : control_component_list})

        elif "save_input" in request.POST and \
                request.POST["save_input"]=="Save input":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "control_id" in request.POST:
                control_id = request.POST["control_id"]
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    config_control_file = sim_para_model.controlfile_set.get(id=control_id)

                    received_input_form = models.ControlInputsForm(request.POST)
                    if received_input_form.is_valid():
                        received_input_data = received_input_form.cleaned_data
                        new_control_input = models.ControlInputs()
                        new_control_input.input_variable_name = \
                                received_input_data["input_variable_name"]
                        if "input_source" in request.POST:
                            new_control_input.input_source = request.POST["input_source"]
                        new_control_input.control_file = config_control_file
                        new_control_input.save()
                        config_control_file.save()
                        sim_para_model.save()

                    try:
                        meter_list = sim_para_model.metercomponents_set.all()
                    except:
                        meter_list = []

                    meter_name_list = []
                    for meter_item in meter_list:
                        meter_name_list.append(meter_item)

                    control_component_list = generating_control_comp_list(\
                                    sim_para_model, config_control_file, 0)

                    return render(request,
                        "config_control_files.html",
                        {'sim_id' : sim_id,
                        'sim_state' : sim_state,
                        'control_id' : control_id,
                        'control_component_list' : control_component_list,
                        'meter_list' : meter_name_list})

        elif "save_output" in request.POST and \
                request.POST["save_output"]=="Save output":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "control_id" in request.POST:
                control_id = request.POST["control_id"]
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    config_control_file = sim_para_model.controlfile_set.get(id=control_id)

                    received_output_form = models.ControlOutputsForm(request.POST)
                    if received_output_form.is_valid():
                        received_output_data = received_output_form.cleaned_data
                        new_control_output = models.ControlOutputs()
                        new_control_output.output_variable_name = \
                                received_output_data["output_variable_name"]
                        new_control_output.output_initial_value = \
                                received_output_data["output_initial_value"]
                        if "output_target" in request.POST:
                            new_control_output.output_target = request.POST["output_target"]
                        new_control_output.control_file = config_control_file
                        new_control_output.save()
                        config_control_file.save()
                        sim_para_model.save()

                    try:
                        meter_list = sim_para_model.metercomponents_set.all()
                    except:
                        meter_list = []

                    meter_name_list = []
                    for meter_item in meter_list:
                        meter_name_list.append(meter_item)

                    control_component_list = generating_control_comp_list(\
                                    sim_para_model, config_control_file, 0)

                    return render(request,
                        "config_control_files.html",
                        {'sim_id' : sim_id,
                        'sim_state' : sim_state,
                        'control_id' : control_id,
                        'control_component_list' : control_component_list,
                        'meter_list' : meter_name_list})

        elif "save_staticvar" in request.POST and \
                request.POST["save_staticvar"]=="Save static variable":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "control_id" in request.POST:
                control_id = request.POST["control_id"]
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    config_control_file = sim_para_model.controlfile_set.get(id=control_id)

                    received_staticvar_form = models.ControlStaticVariableForm(request.POST)
                    if received_staticvar_form.is_valid():
                        received_staticvar_data = received_staticvar_form.cleaned_data
                        new_control_staticvar = models.ControlStaticVariable()
                        new_control_staticvar.static_variable_name = \
                                received_staticvar_data["static_variable_name"]
                        new_control_staticvar.static_initial_value = \
                                received_staticvar_data["static_initial_value"]
                        new_control_staticvar.control_file = config_control_file
                        new_control_staticvar.save()
                        config_control_file.save()
                        sim_para_model.save()

                    try:
                        meter_list = sim_para_model.metercomponents_set.all()
                    except:
                        meter_list = []

                    meter_name_list = []
                    for meter_item in meter_list:
                        meter_name_list.append(meter_item)

                    control_component_list = generating_control_comp_list(\
                                    sim_para_model, config_control_file, 0)

                    return render(request,
                        "config_control_files.html",
                        {'sim_id' : sim_id,
                        'sim_state' : sim_state,
                        'control_id' : control_id,
                        'control_component_list' : control_component_list})

        elif "save_timeevent" in request.POST and \
                request.POST["save_timeevent"]=="Save time event":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "control_id" in request.POST:
                control_id = request.POST["control_id"]
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    config_control_file = sim_para_model.controlfile_set.get(id=control_id)

                    received_timeevent_form = models.ControlTimeEventForm(request.POST)
                    if received_timeevent_form.is_valid():
                        received_timeevent_data = received_timeevent_form.cleaned_data
                        new_control_timeevent = models.ControlTimeEvent()
                        new_control_timeevent.time_event_name = \
                                received_timeevent_data["time_event_name"]
                        new_control_timeevent.initial_time_value = \
                                received_timeevent_data["initial_time_value"]
                        new_control_timeevent.control_file = config_control_file
                        new_control_timeevent.save()
                        config_control_file.save()
                        sim_para_model.save()

                    try:
                        meter_list = sim_para_model.metercomponents_set.all()
                    except:
                        meter_list = []

                    meter_name_list = []
                    for meter_item in meter_list:
                        meter_name_list.append(meter_item)

                    control_component_list = generating_control_comp_list(\
                                    sim_para_model, config_control_file, 0)

                    return render(request,
                        "config_control_files.html",
                        {'sim_id' : sim_id,
                        'sim_state' : sim_state,
                        'control_id' : control_id,
                        'control_component_list' : control_component_list})

        elif "save_varstore" in request.POST and \
                request.POST["save_varstore"]=="Save variable storage":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "control_id" in request.POST:
                control_id = request.POST["control_id"]
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    config_control_file = sim_para_model.controlfile_set.get(id=control_id)

                    received_varstore_form = models.ControlVariableStorageForm(request.POST)
                    if received_varstore_form.is_valid():
                        received_varstore_data = received_varstore_form.cleaned_data
                        new_control_varstore = models.ControlVariableStorage()
                        new_control_varstore.variable_storage_name = \
                                received_varstore_data["variable_storage_name"]
                        new_control_varstore.storage_initial_value = \
                                received_varstore_data["storage_initial_value"]
                        new_control_varstore.storage_status = \
                                received_varstore_data["storage_status"]
                        new_control_varstore.control_file_name = \
                                config_control_file.control_file_name
                        new_control_varstore.sim_case = sim_para_model
                        new_control_varstore.save()
                        config_control_file.save()
                        sim_para_model.save()

                    try:
                        meter_list = sim_para_model.metercomponents_set.all()
                    except:
                        meter_list = []

                    meter_name_list = []
                    for meter_item in meter_list:
                        meter_name_list.append(meter_item)

                    control_component_list = generating_control_comp_list(\
                                    sim_para_model, config_control_file, 0)

                    return render(request,
                        "config_control_files.html",
                        {'sim_id' : sim_id,
                        'sim_state' : sim_state,
                        'control_id' : control_id,
                        'control_component_list' : control_component_list})

        else:
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    ckt_file_list = sim_para_model.circuitschematics_set.all()

                    ckt_ids_for_removal = []
                    for ckt_file_item in ckt_file_list:
                        ckt_ids_for_removal.append("change_ckt_id"+"_"+str(ckt_file_item.id))

                    if ckt_ids_for_removal:
                        for ckt_ids in ckt_ids_for_removal:
                            if ckt_ids in request.POST and request.POST[ckt_ids]=="Remove circuit":
                                for ckt_item_id in ckt_ids_for_removal:
                                    if ckt_item_id in request.POST:
                                        del_ckt_id = int(ckt_item_id.split("_")[-1])
                                        del_ckt = CircuitSchematics.objects.get(id=del_ckt_id)
                                        del_ckt.delete()
                                        sim_para_model.save()

                                simulation_form = []
                                simulation_form.append(sim_para_model)
                                simulation_form.append([])
                                ckt_file_list = sim_para_model.circuitschematics_set.all()
                                ckt_schematic_form = []
                                if ckt_file_list:
                                    for ckt_file_item in ckt_file_list:
                                        ckt_schematic_form.append([ckt_file_item, \
                                                    CircuitSchematicsForm(instance=ckt_file_item)])
                                else:
                                    ckt_schematic_form.append([[], CircuitSchematicsForm()])

                                return render(request,
                                    "edit_circuit_schematic.html",
                                    {'sim_id' : sim_id,
                                    'sim_state' : sim_state,
                                    'ckt_schematic_form' : ckt_schematic_form})

                    try:
                        control_file_list = sim_para_model.controlfile_set.all()
                    except:
                        control_file_list = []

                    control_ids_for_removal = []
                    for control_file_item in control_file_list:
                        control_ids_for_removal.append("change_control_id" + \
                                "_" + str(control_file_item.id))

                    if control_ids_for_removal:
                        for control_ids in control_ids_for_removal:
                            if control_ids in request.POST and \
                                    request.POST[control_ids]=="Remove control":
                                del_control_id = int(control_ids.split("_")[-1])
                                del_control = models.ControlFile.objects.get(\
                                        id=del_control_id)
                                del_control.delete()
                                sim_para_model.save()

                                control_file_list = sim_para_model.controlfile_set.all()
                                control_file_form = []
                                if control_file_list:
                                    for control_file_item in control_file_list:
                                        control_file_form.append([control_file_item, \
                                                    models.ControlFileForm(\
                                                    instance=control_file_item)])
                                else:
                                    control_file_form.append([[], models.ControlFileForm()])

                                ckt_errors = 0
                                ckt_error_list = []

                                return render(request,
                                    "add_control_files.html",
                                    {'sim_id' : sim_id,
                                    'sim_state' : sim_state,
                                    'control_file_form' : control_file_form,
                                    'ckt_errors' : ckt_errors,
                                    'ckt_error_list' : ckt_error_list})

                    components_found, component_objects, \
                            circuit_analysis_components, ckt_error_list = \
                            generate_component_data(sim_para_model)

                    ckts_for_para_update = []
                    for ckt_file_item in ckt_file_list:
                        ckts_for_para_update.append("edit_ckt_para"+"_"+str(ckt_file_item.id))
                    if ckts_for_para_update:
                        for ckt_ids in ckts_for_para_update:
                            if ckt_ids in request.POST and request.POST[ckt_ids]=="View components":
                                ckt_component_list = []
                                recd_ckt_id = int(ckt_ids.split("_")[-1])
                                recd_ckt_item = CircuitSchematics.objects.get(id=recd_ckt_id)
                                if not ckt_error_list:
                                    branch_map = circuit_analysis_components[1]
                                    for comp_items in component_objects.keys():
                                        comp_form_data = []
                                        comp_info = component_objects[comp_items].\
                                                comp_as_a_dict(recd_ckt_item)
                                        comp_model = component_objects[comp_items].\
                                                list_existing_components(recd_ckt_item)
                                        comp_form = component_objects[comp_items].\
                                                comp_as_a_form(recd_ckt_item)
                                        if comp_info:
                                            comp_form_data.append(comp_info)
                                            comp_form_data.append(comp_model)
                                            comp_form_data.append(comp_form)
                                            comp_form_data.append(1)
                                            ckt_component_list.append(comp_form_data)

                                return render(request,
                                    "edit_circuit_parameters.html",
                                    {'sim_id' : sim_id,
                                    'sim_state' : sim_state,
                                    'ckt_component_list' : ckt_component_list,
                                    'ckt_errors' : ckt_error_list})

                    comps_para_update = []
                    for ckt_file_item in ckt_file_list:
                        for comp_items in component_objects.keys():
                            comp_model = component_objects[comp_items].\
                                    list_existing_components(ckt_file_item)
                            if comp_model:
                                comps_para_update.\
                                        append("edit_comp_para"+"_"+str(comp_model.comp_pos_3D))

                    if comps_para_update:
                        for comp_ids in comps_para_update:
                            if comp_ids in request.POST and request.POST[comp_ids]=="Edit parameters":
                                recd_comp_pos3D = comp_ids.split("_")[-1]
                                for ckt_file_item in ckt_file_list:
                                    recd_comp_item = component_objects[recd_comp_pos3D].\
                                                list_existing_components(ckt_file_item)
                                    if recd_comp_item:
                                        break

                                recd_ckt_item = recd_comp_item.comp_ckt
                                ckt_component_list = []
                                if not ckt_error_list:
                                    branch_map = circuit_analysis_components[1]
                                    for comp_items in component_objects.keys():
                                        comp_form_data = []
                                        comp_info = component_objects[comp_items].\
                                                comp_as_a_dict(recd_ckt_item)
                                        comp_model = component_objects[comp_items].\
                                                list_existing_components(recd_ckt_item)
                                        comp_form = component_objects[comp_items].\
                                                comp_as_a_form(recd_ckt_item)
                                        if comp_info:
                                            comp_form_data.append(comp_info)
                                            comp_form_data.append(comp_model)
                                            comp_form_data.append(comp_form)
                                            if comp_model==recd_comp_item:
                                                comp_form_data.append(0)
                                            else:
                                                comp_form_data.append(1)
                                            ckt_component_list.append(comp_form_data)

                                return render(request,
                                    "edit_circuit_parameters.html",
                                    {'sim_id' : sim_id,
                                    'sim_state' : sim_state,
                                    'ckt_component_list' : ckt_component_list,
                                    'ckt_errors' : ckt_error_list})
                    
                    comps_para_submit = []
                    for ckt_file_item in ckt_file_list:
                        for comp_items in component_objects.keys():
                            comp_model = component_objects[comp_items].\
                                    list_existing_components(ckt_file_item)
                            if comp_model:
                                comps_para_submit.\
                                        append("submit_comp_para"+"_"+str(comp_model.comp_pos_3D))

                    if comps_para_submit:
                        for comp_ids in comps_para_submit:
                            if comp_ids in request.POST and request.POST[comp_ids]=="Save parameters":
                                recd_comp_pos3D = comp_ids.split("_")[-1]
                                for ckt_file_item in ckt_file_list:
                                    recd_comp_item = component_objects[recd_comp_pos3D].\
                                                list_existing_components(ckt_file_item)
                                    if recd_comp_item:
                                        form_status = component_objects[recd_comp_pos3D].\
                                                update_form_data(request, recd_comp_item, \
                                                circuit_analysis_components[1])
                                        ckt_file_item.save()
                                        sim_para_model.save()
                                        break

                                recd_ckt_item = recd_comp_item.comp_ckt
                                ckt_component_list = []
                                if not ckt_error_list:
                                    branch_map = circuit_analysis_components[1]
                                    for comp_items in component_objects.keys():
                                        comp_form_data = []
                                        comp_info = component_objects[comp_items].\
                                                comp_as_a_dict(recd_ckt_item)
                                        comp_model = component_objects[comp_items].\
                                                list_existing_components(recd_ckt_item)
                                        comp_form = component_objects[comp_items].\
                                                comp_as_a_form(recd_ckt_item)
                                        if comp_info:
                                            comp_form_data.append(comp_info)
                                            comp_form_data.append(comp_model)
                                            if comp_model==recd_comp_item and form_status:
                                                comp_form_data.append(form_status[0])
                                                comp_form_data.append(0)
                                            else:
                                                comp_form_data.append(comp_form)
                                                comp_form_data.append(1)
                                            ckt_component_list.append(comp_form_data)

                                return render(request,
                                    "edit_circuit_parameters.html",
                                    {'sim_id' : sim_id,
                                    'sim_state' : sim_state,
                                    'ckt_component_list' : ckt_component_list,
                                    'ckt_errors' : ckt_error_list})

                    generate_plot_list = []
                    for ckt_plot_items in sim_para_model.circuitplot_set.all():
                        generate_plot_list.append("generate_plot" + \
                                "_" + str(ckt_plot_items.id))
                    
                    if generate_plot_list:
                        for generate_plot_item in generate_plot_list:
                            if generate_plot_item in request.POST and \
                                    request.POST[generate_plot_item]=="Plot":
                                ckt_plot_id = int(generate_plot_item.split("_")[-1])
                                ckt_plot_item = sim_para_model.circuitplot_set.get(id=ckt_plot_id)
                                y_var_indices = []
                                y_var = []
                                y_var_labels = []
                                y_var_scale = []
                                for waveform_items in ckt_plot_item.circuitwaveforms_set.all():
                                    y_var_labels.append(waveform_items.waveform_legend)
                                    y_var_scale.append(waveform_items.waveform_scale)
                                    for waveform_plots in waveform_items.waveform.all():
                                        y_var_indices.append(waveform_plots.line_pos)
                                        y_var.append([])

                                outputfile_path = os.path.join(os.sep, \
                                    sim_para_model.sim_working_directory, \
                                    sim_para_model.sim_output_file)

                                # Try to read the file.
                                try:
                                    outputfile = open(outputfile_path, "r")
                                except:
                                    pass
                                else:
                                    x_var = []
                                    for line_output in outputfile:
                                        try:
                                            xvalue = line_output.split()[0]
                                        except:
                                            pass
                                        else:
                                            x_var.append(xvalue)
                                            
                                        for c1 in range(len(y_var)):
                                            try:
                                                yvalue = float(line_output.split()[y_var_indices[c1]])
                                            except:
                                                pass
                                            else:
                                                y_var[c1].append(yvalue*y_var_scale[c1])
                                    
                                    for c1 in range(len(x_var)-1, len(x_var)-15, -1):
                                        del x_var[c1]
                                    
                                    for c1 in range(len(y_var)):
                                        for c2 in range(len(y_var[c1])-1, len(x_var)-1, -1):
                                            del y_var[c1][c2]
                                
                                figoutputname = ""
                                for c1 in ckt_plot_item.plot_title.split():
                                    figoutputname += c1
                                figoutputname += ".png"
                                
                                figfile_path = os.path.join(os.sep, \
                                    sim_para_model.sim_working_directory, \
                                    figoutputname)
                                try:
                                    os.remove(figfile_path)
                                except:
                                    pass

                                for c1 in range(len(y_var)):
                                    Matpyplot.plot(x_var, y_var[c1], label=y_var_labels[c1])
                                    Matpyplot.legend()
                                Matpyplot.savefig(figfile_path)
                                Matpyplot.close()
#                                Matpyplot.show()
                                
                                if "sim"+str(sim_para_model.id) in simulation_iteration_collection.keys():
                                    run_state = 1
                                else:
                                    run_state = 0

                                plot_form_list = []
                                for ckt_plots in sim_para_model.circuitplot_set.all():
                                    prev_plot_list = []
                                    for plot_items in ckt_plots.circuitwaveforms_set.all():
                                        prev_plot_list.append([plot_items, 0])
                                    plot_form_list.append([[ckt_plots, 0], prev_plot_list, 0])
                                
                                ckt_error_list = []

                                return render(request,
                                    "output_interface.html",
                                    {'sim_id' : sim_id,
                                    'sim_state' : sim_state,
                                    'run_state' : run_state,
                                    'plot_id' : -1,
                                    'ckt_error_list' : ckt_error_list,
                                    'plot_form_list' : plot_form_list})

                    delete_plot_list = []
                    for ckt_plot_items in sim_para_model.circuitplot_set.all():
                        delete_plot_list.append("delete_plot" + \
                                "_" + str(ckt_plot_items.id))

                    if delete_plot_list:
                        for delete_plot_item in delete_plot_list:
                            if delete_plot_item in request.POST and \
                                    request.POST[delete_plot_item]=="Delete plot":
                                ckt_plot_id = int(delete_plot_item.split("_")[-1])
                                ckt_plot_item = sim_para_model.circuitplot_set.get(id=ckt_plot_id)
                                ckt_plot_item.delete()
                                sim_para_model.save()

                                if "sim"+str(sim_para_model.id) in simulation_iteration_collection.keys():
                                    run_state = 1
                                else:
                                    run_state = 0

                                plot_form_list = []
                                for ckt_plots in sim_para_model.circuitplot_set.all():
                                    prev_plot_list = []
                                    for plot_items in ckt_plots.circuitwaveforms_set.all():
                                        prev_plot_list.append([plot_items, 0])
                                    plot_form_list.append([[ckt_plots, 0], prev_plot_list, 0])
                                
                                ckt_error_list = []

                                return render(request,
                                    "output_interface.html",
                                    {'sim_id' : sim_id,
                                    'sim_state' : sim_state,
                                    'run_state' : run_state,
                                    'plot_id' : -1,
                                    'ckt_error_list' : ckt_error_list,
                                    'plot_form_list' : plot_form_list})

                    config_control_files = []
                    for control_file_item in control_file_list:
                        
                        config_control_files.append("config_control_id" + \
                                "_" + str(control_file_item.id))

                    if config_control_files:
                        for control_ids in config_control_files:
                            if control_ids in request.POST and \
                                    request.POST[control_ids]=="Configure control":
                                config_file_id = int(control_ids.split("_")[-1])
                                config_control_file = models.ControlFile.objects.get(\
                                        id=config_file_id)

                                control_component_list = generating_control_comp_list(\
                                        sim_para_model, config_control_file, 0)

                                return render(request,
                                    "config_control_files.html",
                                    {'sim_id' : sim_id,
                                    'sim_state' : sim_state,
                                    'control_id' : config_file_id,
                                    'control_component_list' : control_component_list,})

                    if "control_id" in request.POST:
                        control_id = int(request.POST["control_id"])
                        control_file_item = sim_para_model.controlfile_set.get(id=control_id)
                        
                        control_inputs_to_delete = []
                        for input_items in control_file_item.controlinputs_set.all():
                            control_inputs_to_delete.append("delete_input" + \
                                    "_" + str(input_items.id))

                        if control_inputs_to_delete:
                            for control_input_ids in control_inputs_to_delete:
                                if control_input_ids in request.POST and \
                                        request.POST[control_input_ids]=="Delete":
                                    input_id = int(control_input_ids.split("_")[-1])
                                    input_item = control_file_item.controlinputs_set.get(\
                                            id=input_id)
                                    input_item.delete()
                                    control_file_item.save()
                                    sim_para_model.save()

                                    control_component_list = generating_control_comp_list(\
                                                sim_para_model, control_file_item, 0)

                                    return render(request,
                                        "config_control_files.html",
                                        {'sim_id' : sim_id,
                                        'sim_state' : sim_state,
                                        'control_id' : control_id,
                                        'control_component_list' : control_component_list,})

                        control_outputs_to_delete = []
                        for output_items in control_file_item.controloutputs_set.all():
                            control_outputs_to_delete.append("delete_output" + \
                                    "_" + str(output_items.id))

                        if control_outputs_to_delete:
                            for control_output_ids in control_outputs_to_delete:
                                if control_output_ids in request.POST and \
                                        request.POST[control_output_ids]=="Delete":
                                    output_id = int(control_output_ids.split("_")[-1])
                                    output_item = control_file_item.controloutputs_set.get(\
                                            id=output_id)
                                    output_item.delete()
                                    control_file_item.save()
                                    sim_para_model.save()

                                    control_component_list = generating_control_comp_list(\
                                                sim_para_model, control_file_item, 0)

                                    return render(request,
                                        "config_control_files.html",
                                        {'sim_id' : sim_id,
                                        'sim_state' : sim_state,
                                        'control_id' : control_id,
                                        'control_component_list' : control_component_list,})

                        control_staticvars_to_delete = []
                        for staticvar_items in control_file_item.controlstaticvariable_set.all():
                            control_staticvars_to_delete.append("delete_staticvar" + \
                                    "_" + str(staticvar_items.id))

                        if control_staticvars_to_delete:
                            for control_staticvar_ids in control_staticvars_to_delete:
                                if control_staticvar_ids in request.POST and \
                                        request.POST[control_staticvar_ids]=="Delete":
                                    staticvar_id = int(control_staticvar_ids.split("_")[-1])
                                    staticvar_item = control_file_item.controlstaticvariable_set.\
                                            get(id=staticvar_id)
                                    staticvar_item.delete()
                                    control_file_item.save()
                                    sim_para_model.save()

                                    control_component_list = generating_control_comp_list(\
                                                sim_para_model, control_file_item, 0)

                                    return render(request,
                                        "config_control_files.html",
                                        {'sim_id' : sim_id,
                                        'sim_state' : sim_state,
                                        'control_id' : control_id,
                                        'control_component_list' : control_component_list,})

                        control_timeevent_to_delete = []
                        for timeevent_items in control_file_item.controltimeevent_set.all():
                            control_timeevent_to_delete.append("delete_timeevent" + \
                                    "_" + str(timeevent_items.id))

                        if control_timeevent_to_delete:
                            for control_timeevent_ids in control_timeevent_to_delete:
                                if control_timeevent_ids in request.POST and \
                                        request.POST[control_timeevent_ids]=="Delete":
                                    timeevent_id = int(control_timeevent_ids.split("_")[-1])
                                    timeevent_item = control_file_item.controltimeevent_set.\
                                            get(id=timeevent_id)
                                    timeevent_item.delete()
                                    control_file_item.save()
                                    sim_para_model.save()

                                    control_component_list = generating_control_comp_list(\
                                                sim_para_model, control_file_item, 0)

                                    return render(request,
                                        "config_control_files.html",
                                        {'sim_id' : sim_id,
                                        'sim_state' : sim_state,
                                        'control_id' : control_id,
                                        'control_component_list' : control_component_list,})

                        control_varstore_to_delete = []
                        for varstore_items in sim_para_model.controlvariablestorage_set.\
                                    filter(control_file_name=control_file_item.control_file_name):
                            control_varstore_to_delete.append("delete_varstore" + \
                                    "_" + str(varstore_items.id))

                        if control_varstore_to_delete:
                            for control_varstore_ids in control_varstore_to_delete:
                                if control_varstore_ids in request.POST and \
                                        request.POST[control_varstore_ids]=="Delete":
                                    varstore_id = int(control_varstore_ids.split("_")[-1])
                                    varstore_item = sim_para_model.controlvariablestorage_set.\
                                            get(id=varstore_id)
                                    varstore_item.delete()
                                    sim_para_model.save()

                                    control_component_list = generating_control_comp_list(\
                                                sim_para_model, control_file_item, 0)

                                    return render(request,
                                        "config_control_files.html",
                                        {'sim_id' : sim_id,
                                        'sim_state' : sim_state,
                                        'control_id' : control_id,
                                        'control_component_list' : control_component_list,})

            choosing_simulations = []
            sim_list = SimulationCase.objects.all()
            for sim_item in sim_list:
                choosing_simulations.append("choose_sim"+"_"+str(sim_item.id))
            
            for sim_item in choosing_simulations:
                print(sim_item)
                if sim_item in request.POST and request.POST[sim_item]=="Load simulation":
                    sim_id = int(sim_item.split("_")[-1])

                    sim_state, simulation_form = load_simulation_parameters(request, sim_id)
                    return render(request,
                            "edit_simulation.html",
                            {'sim_id' : sim_id,
                            'sim_state' : sim_state,
                            'simulation_form' : simulation_form})

    if sim_state==1:
        return render(request,
                "edit_simulation.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'simulation_form' : simulation_form})
    
    elif sim_state==2:
        return render(request,
                "edit_circuit_schematic.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'ckt_schematic_form' : ckt_schematic_form})
    
    else:
        return render(request,
                "edit_simulation.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'simulation_form' : simulation_form})



def simulation_library(request):
    simulation_collection = SimulationCase.objects.all()
    return render(request,
                "list_simulation.html",
                {'simulation_collection' : simulation_collection, })
