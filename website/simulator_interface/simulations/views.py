from django.shortcuts import render
from django.http import HttpResponse
from models import SimulationCase, SimulationCaseForm, CircuitSchematics, CircuitSchematicsForm

# Create your views here.

def index(request):
    return render(request, "index.html")


def new_simulation(request):
    if not request.method == "POST":
        simulation_form = []
        simulation_form.append(SimulationCaseForm())
        ckt_schematic_form = []
        ckt_schematic_form.append([[], CircuitSchematicsForm()])

    else:
        if "submit" in request.POST and request.POST["submit"]=="Submit":
            sim_para_received = SimulationCaseForm(request.POST)
            if sim_para_received.is_valid():
                sim_parameters = sim_para_received.cleaned_data
                sim_para_model = SimulationCase()
                sim_para_model.sim_title = sim_parameters["sim_title"]
                sim_para_model.sim_descrip = sim_parameters["sim_descrip"]
                sim_para_model.sim_time_limit = sim_parameters["sim_time_limit"]
                sim_para_model.sim_time_data = sim_parameters["sim_time_data"]
                sim_para_model.sim_output_file = sim_parameters["sim_output_file"]
                sim_para_model.sim_output_slice = sim_parameters["sim_output_slice"]
                sim_para_model.sim_div_number = sim_parameters["sim_div_number"]
                sim_para_model.sim_working_directory = sim_parameters["sim_working_directory"]
                
                ckt_file = CircuitSchematics()
                if u'ckt_file_path' in request.POST:
                    ckt_file.ckt_file_name = request.POST.getlist(u'ckt_file_path')[0]
                if u'ckt_file_descrip' in request.POST:
                    ckt_file.ckt_file_descrip = request.POST.getlist(u'ckt_file_descrip')[0]
                ckt_file.save()
                sim_para_model.sim_circuit_files = ckt_file
                sim_para_model.save()
            
                simulation_form = []
                simulation_form.append(SimulationCaseForm(instance=sim_para_model))
                ckt_schematic_form = []
                ckt_schematic_form.append([ckt_file, CircuitSchematicsForm(instance=ckt_file)])

    return render(request,
          "new_simulation.html",
          {'simulation_form' : simulation_form,
           'ckt_schematic_form' : ckt_schematic_form})
