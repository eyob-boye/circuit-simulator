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
    if sim_para_received.is_valid():
        sim_parameters = sim_para_received.cleaned_data
        if "sim_id" in request.POST:
            sim_id = int(request.POST["sim_id"])
            if sim_id>0:
                sim_para_model = SimulationCase.objects.get(id=sim_id)
            else:
                sim_para_model = SimulationCase()
        sim_para_model.sim_title = sim_parameters["sim_title"]
        sim_para_model.sim_descrip = sim_parameters["sim_descrip"]
        sim_para_model.sim_time_limit = sim_parameters["sim_time_limit"]
        sim_para_model.sim_time_data = sim_parameters["sim_time_data"]
        sim_para_model.sim_output_file = sim_parameters["sim_output_file"]
        sim_para_model.sim_output_slice = sim_parameters["sim_output_slice"]
        sim_para_model.sim_div_number = sim_parameters["sim_div_number"]
        sim_para_model.sim_working_directory = sim_parameters["sim_working_directory"]
        sim_para_model.save()

        simulation_form = []
        simulation_form.append([])
        simulation_form.append([SimulationCaseForm(instance=sim_para_model)])
    else:
        pass


def index(request):
    return render(request, "index.html")


def new_simulation(request):
    if not request.method == "POST":
        sim_id = -1
        simulation_form = []
        simulation_form.append(SimulationCaseForm())
        ckt_schematic_form = []
        ckt_schematic_form.append([[], CircuitSchematicsForm()])

    else:
        print(request.POST)
        print(repr(request.POST["sim_working_directory"]))
        print
        print
        if "add_ckt_schematic" in request.POST and request.POST["add_ckt_schematic"]=="Add circuit file":
            sim_para_received = SimulationCaseForm(request.POST)
            if sim_para_received.is_valid():
                sim_parameters = sim_para_received.cleaned_data
                if "sim_id" in request.POST:
                    sim_id = int(request.POST["sim_id"])
                    if sim_id>0:
                        sim_para_model = SimulationCase.objects.get(id=sim_id)
                    else:
                        sim_para_model = SimulationCase()
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

                if u'ckt_file_path' in request.POST:
                    ckt_file = CircuitSchematics()
                    ckt_file.ckt_file_name = request.POST.getlist(u'ckt_file_path')[0]
                    if u'ckt_file_descrip' in request.POST:
                        ckt_file.ckt_file_descrip = request.POST.getlist(u'ckt_file_descrip')[0]
                    ckt_file.ckt_sim_case = sim_para_model
                    ckt_file.save()


                simulation_form = []
                simulation_form.append(SimulationCaseForm(instance=sim_para_model))
                ckt_file_list = sim_para_model.circuitschematics_set.all()
                ckt_schematic_form = []
                for ckt_file_item in ckt_file_list:
                    ckt_schematic_form.append([ckt_file_item, CircuitSchematicsForm(instance=ckt_file_item)])
                ckt_schematic_form.append([[], CircuitSchematicsForm()])

        elif "submit" in request.POST and request.POST["submit"]=="Submit":
            sim_para_received = SimulationCaseForm(request.POST)
            if sim_para_received.is_valid():
                sim_parameters = sim_para_received.cleaned_data
                if "sim_id" in request.POST:
                    sim_id = int(request.POST["sim_id"])
                    if sim_id>0:
                        sim_para_model = SimulationCase.objects.get(id=sim_id)
                    else:
                        sim_para_model = SimulationCase()
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

                if u'ckt_file_path' in request.POST:
                    ckt_file = CircuitSchematics()
                    ckt_file.ckt_file_name = request.POST.getlist(u'ckt_file_path')[0]
                    if u'ckt_file_descrip' in request.POST:
                        ckt_file.ckt_file_descrip = request.POST.getlist(u'ckt_file_descrip')[0]
                    ckt_file.ckt_sim_case = sim_para_model
                    ckt_file.save()


                simulation_form = []
                simulation_form.append(SimulationCaseForm(instance=sim_para_model))
                ckt_file_list = sim_para_model.circuitschematics_set.all()
                ckt_schematic_form = []
                for ckt_file_item in ckt_file_list:
                    ckt_schematic_form.append([ckt_file_item, CircuitSchematicsForm(instance=ckt_file_item)])

                if "\\" in sim_para_model.sim_working_directory:
                    os_type = "win"
                elif "/" in sim_para_model.sim_working_directory:
                    os_type = "unix"

                f = open(os.path.join(os.sep, \
                                      sim_para_model.sim_working_directory, \
                                      ckt_file_list[0].ckt_file_name), "r")
                for line in f:
                    print(line)
            
            else:
                sim_id = -1
                simulation_form = []
                simulation_form.append(SimulationCaseForm(request.POST))
                ckt_schematic_form = []
                ckt_schematic_form.append([[], CircuitSchematicsForm()])
                

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
                                sim_para_received = SimulationCaseForm(request.POST)
                                if sim_para_received.is_valid():
                                    sim_parameters = sim_para_received.cleaned_data
                                    if "sim_id" in request.POST:
                                        sim_id = int(request.POST["sim_id"])
                                        if sim_id>0:
                                            sim_para_model = SimulationCase.objects.get(id=sim_id)
                                        else:
                                            sim_para_model = SimulationCase()
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
                                    simulation_form.append(SimulationCaseForm(instance=sim_para_model))
                                    ckt_file_list = sim_para_model.circuitschematics_set.all()
                                    ckt_schematic_form = []
                                    for ckt_file_item in ckt_file_list:
                                        ckt_schematic_form.append([ckt_file_item, CircuitSchematicsForm(instance=ckt_file_item)])

    return render(request,
            "new_simulation.html",
            {'sim_id' : sim_id, 
            'simulation_form' : simulation_form,
            'ckt_schematic_form' : ckt_schematic_form})
