from django.db import models
from django import forms
from django.forms import ModelForm

# Create your models here.


class SimulationCase(models.Model):
    """
    Contains the overall simulation parameters.
    """
    sim_title = models.CharField(max_length=100)
    sim_descrip = models.TextField(blank=True, null=True)
    sim_time_limit = models.FloatField(default=1.0)
    sim_time_step = models.FloatField(default=1.0e-6)
    sim_time_data = models.FloatField(default=10.0e-6)
    sim_output_file = models.CharField(max_length=30, default="ckt_output.dat")
    sim_output_slice = models.CharField(max_length=3,
                                        choices=(("Yes", "Yes"),
                                                 ("No", "No")),
                                        default="No")
    sim_div_number = models.IntegerField(default=1)

    def __unicode__(self):
        return self.sim_title


class SimulationCaseForm(ModelForm):
    class Meta:
        model = SimulationCase
        fields = ('sim_title',
                  'sim_descrip',
                  'sim_time_limit',
                  'sim_time_step',
                  'sim_time_data',
                  'sim_output_file',
                  'sim_output_slice',
                  'sim_div_number')
        widgets = {
            'sim_title': forms.TextInput(attrs={'size': 80}),
            'sim_descrip': forms.Textarea(attrs={'rows': 15, 'cols': 80}),
            }
