from django.shortcuts import render
from django.http import HttpResponse
from django.forms.models import model_to_dict
from models import SimulationCase, SimulationCaseForm, CircuitSchematics, CircuitSchematicsForm
import os
import network_reader as NwRdr

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
                    ckt_schematic_form.append([ckt_file_item, CircuitSchematicsForm(instance=ckt_file_item)])
            # If no circuit files, create a blank circuit add form
            else:
                ckt_schematic_form.append([[], CircuitSchematicsForm()])
    
    return [sim_id, sim_state, simulation_form, ckt_schematic_form]


def check_circuit_errors(sim_para_model):
    ckt_file_list = sim_para_model.circuitschematics_set.all()
    ckt_schematic_form = []
    ckt_errors = 0
    # List the existing circuits.
    for ckt_file_item in ckt_file_list:
        ckt_item_dict = model_to_dict(ckt_file_item)
        ckt_item_form = CircuitSchematicsForm(ckt_item_dict, instance=ckt_file_item)
        ckt_full_path = os.path.join(os.sep, sim_para_model.sim_working_directory, ckt_file_item.ckt_file_name)

        # Try to read the file.
        try:
            check_ckt_file = open(ckt_full_path, "r")
        # If can't be read, it means file doesn't exist in the working directory.
        except:

            if ckt_item_form.is_valid():
                ckt_item_form.add_error('ckt_file_descrip', 'Circuit spreadsheet could not be read. \
                                                Make sure it is in same directory as working directory above')
                ckt_errors = 1

            ckt_schematic_form.append([ckt_file_item, ckt_item_form])

        # If it can be read, it could be a genuine file in which case save it.
        # Or else, it may not be a .csv file, in which case raise an error.
        else:
            if len(ckt_file_item.ckt_file_name.split("."))>1 and ckt_file_item.ckt_file_name.split(".")[-1]=="csv":
                ckt_schematic_form.append([ckt_file_item, CircuitSchematicsForm(instance=ckt_file_item)])
            else:
                ckt_item_form.add_error('ckt_file_descrip', 'Circuit schematic must be a .csv file.')
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
        ckt_full_path = os.path.join(os.sep, sim_para_model.sim_working_directory, received_ckt_file_name)
        # Try to read the file.
        try:
            check_ckt_file = open(ckt_full_path, "r")
        # If can't be read, it means file doesn't exist in the working directory.
        except:
            ckt_form = CircuitSchematicsForm(request.POST)
            if ckt_form.is_valid():
                ckt_form.add_error('ckt_file_descrip', 'Circuit spreadsheet could not be read. \
                                                Make sure it is in same directory as working directory above')
            ckt_schematic_form.append([[], ckt_form])

        # If it can be read, it could be a genuine file in which case save it.
        # Or else, it may not be a .csv file, in which case raise an error.
        else:
            if len(received_ckt_file_name.split("."))>1 and received_ckt_file_name.split(".")[-1]=="csv":
                ckt_file.ckt_file_name = received_ckt_file_name
                if u'ckt_file_descrip' in request.POST:
                    ckt_file.ckt_file_descrip = request.POST.getlist(u'ckt_file_descrip')[0]
                ckt_file.ckt_sim_case = sim_para_model
                ckt_file.save()
                sim_para_model.save()
                ckt_schematic_form.append([ckt_file, CircuitSchematicsForm(instance=ckt_file)])
            else:
                ckt_form = CircuitSchematicsForm(request.POST)
                ckt_form.add_error('ckt_file_descrip', 'Circuit schematic must be a .csv file.')
                ckt_schematic_form.append([[], ckt_form])

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



def process_circuit_schematics(sim_para_model):
    ckt_file_list = sim_para_model.circuitschematics_set.all()
    nw_input = []
    conn_ckt_mat = []
    if ckt_file_list:
        for ckt_file_item in ckt_file_list:
            nw_input.append(ckt_file_item.ckt_file_name.split(".csv")[0])
            full_file_path = os.path.join(os.sep, sim_para_model.sim_working_directory, ckt_file_item.ckt_file_name)
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
        node_list, branch_map, node_branch_errors = NwRdr.determine_nodes_branches(conn_ckt_mat, nw_input)
        if node_branch_errors:
            ckt_error_list.extend(node_branch_errors)
        else:
            circuit_analysis_components.append(node_list)
            circuit_analysis_components.append(branch_map)

    # Segregating the elements in the circuit as those
    # having source voltages, meters and those that can be
    # controlled.
    if not ckt_error_list:
        bundled_list_of_elements = NwRdr.classify_components(components_found, component_objects)
        circuit_analysis_components.append(bundled_list_of_elements)
        source_list = bundled_list_of_elements[0]
#        meter_list = bundled_list_of_elements[1]
#        controlled_elements = bundled_list_of_elements[2]

    if not ckt_error_list:
        # Make a list of all the loops in the circuit.
        loop_list, loop_branches = NwRdr.determine_loops(conn_ckt_mat, node_list, branch_map)
        circuit_analysis_components.append(loop_list)
        circuit_analysis_components.append(loop_branches)
        print(loop_list)
        print
        print
        print(loop_branches)
        print
        print

    if not ckt_error_list:
        # Convert the above list of loops as segments of branches
        # while branch_params lists out the branch segments with
        # their parameters as the last element
        system_loops, branch_params = NwRdr.update_branches_loops(loop_branches, source_list)
        circuit_analysis_components.append(system_loops)
        circuit_analysis_components.append(branch_params)

    if not ckt_error_list:
        # Make lists of all the nodes connected together by 
        # empty branches.
        shortnode_list, shortbranch_list = NwRdr.delete_empty_branches(system_loops, branch_params, node_list, component_objects)
        circuit_analysis_components.append(shortnode_list)
        circuit_analysis_components.append(shortbranch_list)

    if not ckt_error_list:
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

        circuit_analysis_components.append(components_in_branch)

        for c1 in components_in_branch:
            print(c1)

    if not ckt_error_list:
        # Branch currents for nodal analysis
        branch_currents = []
        for c1 in range(len(branch_params)):
            branch_currents.append(0.0)
        circuit_analysis_components.append(branch_currents)

    for c1 in ckt_error_list:
        print(c1)
    print
    print

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

        if "save_ckt_schematic" in request.POST and request.POST["save_ckt_schematic"]=="Save circuit file":
            sim_id, sim_state, ckt_schematic_form, ckt_errors = save_circuit_schematic(request)
            if sim_state==1:
                sim_state = 2
            return render(request,
                "edit_circuit_schematic.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'ckt_schematic_form' : ckt_schematic_form,
                'ckt_errors' : ckt_errors})
            

        elif "add_ckt_schematic" in request.POST and request.POST["add_ckt_schematic"]=="Add circuit file":
            sim_id, sim_state, ckt_schematic_form, ckt_errors = add_circuit_schematic(request)
            if sim_state==1:
                sim_state = 2
            return render(request,
                "edit_circuit_schematic.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'ckt_schematic_form' : ckt_schematic_form,
                'ckt_errors' : ckt_errors})

        elif "save_sim_param" in request.POST and request.POST["save_sim_param"]=="Save Simulation Parameters":
            sim_id, sim_state, simulation_form, ckt_schematic_form = save_simulation_parameters(request)
            if sim_state==1:
                sim_state = 2
            return render(request,
                    "edit_simulation.html",
                    {'sim_id' : sim_id,
                    'sim_state' : sim_state,
                    'simulation_form' : simulation_form})

        elif "edit_sim_param" in request.POST and request.POST["edit_sim_param"]=="Edit Simulation Parameters":
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
                            full_file_path = os.path.join(os.sep, sim_para_model.sim_working_directory, ckt_file_item.ckt_file_name)
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

        elif "process_ckt_schematics" in request.POST and request.POST["process_ckt_schematics"]=="Process circuit schematics":
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

            ckt_errors = 0
            ckt_file_list = sim_para_model.circuitschematics_set.all()
            ckt_schematic_form = []
            if ckt_file_list:
                for ckt_file_item in ckt_file_list:
                    ckt_schematic_form.append([ckt_file_item, CircuitSchematicsForm(instance=ckt_file_item)])

            components_found, component_objects, \
                    circuit_analysis_components, ckt_error_list = \
                    process_circuit_schematics(sim_para_model)

            return render(request,
                "edit_circuit_schematic.html",
                {'sim_id' : sim_id,
                'sim_state' : sim_state,
                'ckt_schematic_form' : ckt_schematic_form,
                'ckt_errors' : ckt_errors})

        else:
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
                                if "sim_state" in request.POST:
                                    sim_state = int(request.POST["sim_state"])

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
                                        ckt_schematic_form.append([ckt_file_item, CircuitSchematicsForm(instance=ckt_file_item)])
                                else:
                                    ckt_schematic_form.append([[], CircuitSchematicsForm()])

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
