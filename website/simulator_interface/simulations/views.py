from django.shortcuts import render
from django.http import HttpResponse
from models import SimulationCase, SimulationCaseForm

# Create your views here.

def index(request):
    return render(request, "index.html")


def new_simulation(request):
    simulation_form = SimulationCaseForm()
    return render(request,
                  "new_simulation.html",
                  {'simulation_form' : simulation_form})
