from django.shortcuts import render
from django.http import HttpResponse
from models import SimulationCase, SimulationCaseForm, CircuitSchematics, CircuitSchematicsForm
import os

# Create your views here.


def extract_simulation_case(request):
    """
    This function extracts the form data from a SimulationCase form
    and saves it in a SimulationCase model object.
    """
    sim_para_received = SimulationCaseForm(request.POST)
    if "sim_id" in request.POST:
        sim_id = int(request.POST["sim_id"])
        if sim_id>0:
            sim_para_model = SimulationCase.objects.get(id=sim_id)
        else:
            sim_para_model = SimulationCase()
    if sim_para_received.is_valid():
        sim_parameters = sim_para_received.cleaned_data
#        if "sim_id" in request.POST:
#            sim_id = int(request.POST["sim_id"])
#            if sim_id>0:
#                sim_para_model = SimulationCase.objects.get(id=sim_id)
#            else:
#                sim_para_model = SimulationCase()
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

        simulation_form = []
        simulation_form.append(sim_para_model)
        simulation_form.append([])

    else:
        simulation_form = []
        simulation_form.append([])
        simulation_form.append(sim_para_received)

    return [sim_id, simulation_form]


def index(request):
    return render(request, "index.html")


def new_simulation(request):
    if not request.method == "POST":
        sim_id = -1
        simulation_form = []
        simulation_form.append([])
        simulation_form.append(SimulationCaseForm())
        ckt_schematic_form = []

    else:
        print(request.POST)
        print
        print
        if "save_ckt_schematic" in request.POST and request.POST["save_ckt_schematic"]=="Save circuit file":
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)

            simulation_form = []
            simulation_form.append(sim_para_model)
            simulation_form.append([])
            ckt_file_list = sim_para_model.circuitschematics_set.all()
            ckt_schematic_form = []
            for ckt_file_item in ckt_file_list:
                ckt_schematic_form.append([ckt_file_item, CircuitSchematicsForm(instance=ckt_file_item)])

            if u'ckt_file_path' in request.POST:
                ckt_file = CircuitSchematics()
                received_ckt_file_name = request.POST.getlist(u'ckt_file_path')[0]
                ckt_full_path = os.path.join(os.sep, sim_para_model.sim_working_directory, received_ckt_file_name)
                try:
                    check_ckt_file = open(ckt_full_path, "r")
                except:
                    ckt_form = CircuitSchematicsForm(request.POST)
                    if ckt_form.is_valid():
                        ckt_form.add_error('ckt_file_descrip', 'Circuit spreadsheet could not be read. \
                                                        Make sure in same directory as working directory above')
                    ckt_schematic_form.append([[], ckt_form])
                    print(ckt_form)
                    print("Testing")
                    
                else:
                    ckt_file.ckt_file_name = request.POST.getlist(u'ckt_file_path')[0]
                    if u'ckt_file_descrip' in request.POST:
                        ckt_file.ckt_file_descrip = request.POST.getlist(u'ckt_file_descrip')[0]
                    ckt_file.ckt_sim_case = sim_para_model
                    ckt_file.save()
                    sim_para_model.save()
                    ckt_schematic_form.append([ckt_file, CircuitSchematicsForm(instance=ckt_file)])


        elif "add_ckt_schematic" in request.POST and request.POST["add_ckt_schematic"]=="Add circuit file":
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)

            simulation_form = []
            simulation_form.append(sim_para_model)
            simulation_form.append([])
            ckt_file_list = sim_para_model.circuitschematics_set.all()
            ckt_schematic_form = []
            for ckt_file_item in ckt_file_list:
                ckt_schematic_form.append([ckt_file_item, CircuitSchematicsForm(instance=ckt_file_item)])
            ckt_schematic_form.append([[], CircuitSchematicsForm()])

        elif "save_sim_param" in request.POST and request.POST["save_sim_param"]=="Save Simulation Parameters":
            sim_id, simulation_form = extract_simulation_case(request)
            print(simulation_form[1])
            ckt_schematic_form = []
            if sim_id>0:
                sim_para_model = SimulationCase.objects.get(id=sim_id)
                if not simulation_form[1]:
                    ckt_file_list = sim_para_model.circuitschematics_set.all()
                    if ckt_file_list:
                        for ckt_file_item in ckt_file_list:
                            ckt_schematic_form.append([ckt_file_item, CircuitSchematicsForm(instance=ckt_file_item)])
                    else:
                        ckt_schematic_form.append([[], CircuitSchematicsForm()])

        elif "edit_sim_param" in request.POST and request.POST["edit_sim_param"]=="Edit Simulation Parameters":
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                else:
                    sim_para_model = SimulationCase()

            simulation_form = []
            simulation_form.append([])
            simulation_form.append(SimulationCaseForm(instance=sim_para_model))
            ckt_schematic_form = []

        else:
            if "sim_id" in request.POST:
                sim_id = int(request.POST["sim_id"])
                print(sim_id)
                if sim_id>0:
                    sim_para_model = SimulationCase.objects.get(id=sim_id)
                    ckt_file_list = sim_para_model.circuitschematics_set.all()
                    list_of_ckt_ids = []
                    for ckt_file_item in ckt_file_list:
                        list_of_ckt_ids.append("change_ckt_id"+"_"+str(ckt_file_item.id))                
                    print("Check list")
                    print(list_of_ckt_ids)
                    if list_of_ckt_ids:
                        for ckt_ids in list_of_ckt_ids:
                            if ckt_ids in request.POST and request.POST[ckt_ids]=="Remove circuit":
                                for ckt_item_id in list_of_ckt_ids:
                                    if ckt_item_id in request.POST:
                                        print("Check delete")
                                        print(ckt_item_id)
                                        del_ckt_id = int(ckt_item_id.split("_")[-1])
                                        print(del_ckt_id)
                                        del_ckt = CircuitSchematics.objects.get(id=del_ckt_id)
                                        print(del_ckt.ckt_file_name)
                                        del_ckt.delete()
                                        sim_para_model.save()

                                simulation_form = []
                                simulation_form.append(sim_para_model)
                                simulation_form.append([])
                                ckt_file_list = sim_para_model.circuitschematics_set.all()
                                ckt_schematic_form = []
                                for ckt_file_item in ckt_file_list:
                                    ckt_schematic_form.append([ckt_file_item, CircuitSchematicsForm(instance=ckt_file_item)])

    return render(request,
            "new_simulation.html",
            {'sim_id' : sim_id, 
            'simulation_form' : simulation_form,
            'ckt_schematic_form' : ckt_schematic_form})
