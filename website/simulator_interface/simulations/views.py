from django.shortcuts import render
from django.http import HttpResponse
from django.forms.models import model_to_dict
from models import SimulationCase, SimulationCaseForm, CircuitSchematics, CircuitSchematicsForm, CircuitComponents
import models
import os
import network_reader as NwRdr
import circuit_elements as CktElem
import matrix
import threading


def prepare_simulation_objects(sim_components, ckt_topo, conn_ckt_mat):
    
    synthesized_ckt_comps = {}
    
    components_found = sim_components[0]
    component_objects = sim_components[1]

    node_list = ckt_topo[0]
    branch_map = ckt_topo[1]

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
#
#    import __control
#
#    plotted_variable_list = []
#    for c1 in control_file_variablestorage.keys():
#        if control_file_variablestorage[c1][1]=="yes":
#            plotted_variable_list.append(c1)

    return





def simulation_iterations(sim_id, time):
    sim_para_model = SimulationCase.objects.get(id=int(sim_id))
    f_path = os.path.join(os.sep, \
                sim_para_model.sim_working_directory, \
                "check_threading.txt")
    f = open(f_path, "w")
    t = float(time)
    t_step = 0.1
    t_start = t + t_step
    
    check_path = os.path.join(os.sep, \
                        sim_para_model.sim_working_directory, \
                        "sim_log_run_state.txt")

    while t>0.0:
        try:
            check_status = open(check_path, "r")
        except:
            pass
        else:
            for line in check_status:
                if line=="Start":
                    t = 0.0
                if line=="Run":
                    t += 0.00001
                    if t>=t_start:
                        f.write(str(t))
                        f.write("\n")
                        t_start += t_step
                if line=="Stop":
                    t = 0.0
    check_status.close()
    f.close()
    return

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
                
                f_path = os.path.join(os.sep, \
                                    sim_para_model.sim_working_directory, \
                                    "sim_log_run_state.txt")
                print("c1")
                try:
                    f = open(f_path, "r")
                except:
                    f = open(f_path, "w")
                    f.write("Start")
                    f.close()
                    sim_status = "Start"
                    print("c2")
                else:
                    print("c3")
                    for line in f:
                        if not line=="Run":
                            f.close()
                            f = open(f_path, "w")
                            f.write("Start")
                            f.close()
                        else:
                            sim_status = "Run"

                if sim_status=="Start":
                    simulator_loop = threading.Thread(target=simulation_iterations, \
                            kwargs={'sim_id':sim_id, \
                                    'time':0.0, })
                    simulator_loop.start()
                
                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'ckt_error_list' : ckt_error_list})

        elif "run_simulation" in request.POST and request.POST["run_simulation"]=="Run":
            print("Here")
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)

                    f_path = os.path.join(os.sep, \
                                        sim_para_model.sim_working_directory, \
                                        "sim_log_run_state.txt")
                    f = open(f_path, "w")
                    f.write("Run")
                    f.close()
                
                ckt_error_list = []

                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'ckt_error_list' : ckt_error_list})

        elif "stop_simulation" in request.POST and request.POST["stop_simulation"]=="Stop":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    sim_para_model.sim_run_state = 2
                    sim_para_model.save()

                    f_path = os.path.join(os.sep, \
                                        sim_para_model.sim_working_directory, \
                                        "sim_log_run_state.txt")
                    f = open(f_path, "w")
                    f.write("Stop")
                    f.close()
                
                ckt_error_list = []

                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'ckt_error_list' : ckt_error_list})

        elif "pause_simulation" in request.POST and request.POST["pause_simulation"]=="Pause":
            if "sim_state" in request.POST:
                sim_state = int(request.POST["sim_state"])
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    sim_para_model.sim_run_state = 3
                    sim_para_model.save()
                    f_path = os.path.join(os.sep, \
                                        sim_para_model.sim_working_directory, \
                                        "sim_log_run_state.txt")
                    f = open(f_path, "w")
                    f.write("Pause")
                    f.close()
                    
                ckt_error_list = []

                return render(request,
                    "output_interface.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'ckt_error_list' : ckt_error_list})

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
