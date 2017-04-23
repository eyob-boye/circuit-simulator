from django.shortcuts import render
from django.http import HttpResponse
from models import SimulationCase, SimulationCaseForm

# Create your views here.

def index(request):
    return render(request, "index.html")


def new_simulation(request):
    if not request.method == "POST":
        simulation_form = SimulationCaseForm()
    else:
        if "submit" in request.POST and request.POST["submit"]=="Submit":
            print(request.POST)
            print
            print
            sim_para_received = SimulationCaseForm(request.POST)
            if sim_para_received.is_valid():
                sim_parameters = sim_para_received.cleaned_data
                print(sim_parameters)
                print
                print
                sim_para_model = SimulationCase()
                sim_para_model.sim_title = sim_parameters["sim_title"]
                sim_para_model.sim_descrip = sim_parameters["sim_descrip"]
                sim_para_model.sim_time_limit = sim_parameters["sim_time_limit"]
                sim_para_model.sim_time_data = sim_parameters["sim_time_data"]
                sim_para_model.sim_output_file = sim_parameters["sim_output_file"]
                sim_para_model.sim_output_slice = sim_parameters["sim_output_slice"]
                sim_para_model.sim_div_number = sim_parameters["sim_div_number"]
                
                
            simulation_form = SimulationCaseForm(instance=sim_para_model)
    return render(request,
                  "new_simulation.html",
                  {'simulation_form' : simulation_form})
