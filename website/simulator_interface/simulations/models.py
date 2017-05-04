from django.db import models
from django import forms
from django.forms import ModelForm
import os

# Create your models here.

class SimulationCase(models.Model):
    """
    Contains the overall simulation parameters.
    """
    sim_title = models.CharField(max_length=100, default="Test case", verbose_name="Simulation Title")
    sim_descrip = models.TextField(blank=True, null=True, verbose_name="Simulation description")
    sim_time_limit = models.FloatField(default=1.0, verbose_name="Time duration")
    sim_time_step = models.FloatField(default=1.0e-6, verbose_name="Integration time step")
    sim_time_data = models.FloatField(default=10.0e-6, verbose_name="Time step of data storage")
    sim_output_file = models.CharField(max_length=30, default="ckt_output.dat", verbose_name="Output data file")
    sim_output_slice = models.CharField(max_length=3,
                                        choices=(("Yes", "Yes"),
                                                 ("No", "No")),
                                        default="No",
                                        verbose_name="Slice the output file?")
    sim_div_number = models.IntegerField(default=1,
                                        verbose_name="Number of slices")
    sim_working_directory = models.CharField(max_length=300, verbose_name="Directory with circuit files")

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
                  'sim_div_number', 
                  'sim_working_directory')
        widgets = {
            'sim_title': forms.TextInput(attrs={'size': 80}),
            'sim_descrip': forms.Textarea(attrs={'rows': 15, 'cols': 80}),
            'sim_output_file': forms.TextInput(attrs={'size': 50}),
            'sim_working_directory': forms.TextInput(attrs={'size': 80}),
            }

    def clean_sim_time_limit(self):
        sim_time_limit = float(self.cleaned_data["sim_time_limit"])
        if sim_time_limit<=0.0:
            raise forms.ValidationError("Time limit must be greater than 0.0")
        
        return sim_time_limit


    def clean_sim_working_directory(self):
        webdir = self.cleaned_data["sim_working_directory"]
        try:
            test_file = open(os.path.join(os.sep, \
                            webdir, \
                            "testfile"), "w")
        except:
            raise forms.ValidationError("Can't write into this directory!")
        
        return webdir

    def clean(self):
        cleaned_data = super(SimulationCaseForm, self).clean()
        sim_time_step = float(cleaned_data.get('sim_time_step'))
        sim_time_data = float(cleaned_data.get('sim_time_data'))
        if sim_time_data<sim_time_step:
            self.add_error('sim_time_data', 'This number must be greater than or equal to Integration Time Step (previous)')
        
        sim_output_slice = cleaned_data.get('sim_output_slice')
        sim_div_number = cleaned_data.get('sim_div_number')
        if sim_output_slice=="Yes" and sim_div_number<2:
            self.add_error('sim_div_number', 'If the output file slicing option is chosen, the Number of Slices must be at least 2')


class CircuitSchematics(models.Model):
    ckt_file_path = models.FileField(max_length=300, upload_to='', blank=True)
    ckt_file_descrip = models.CharField(max_length=100, default="Sample circuit", verbose_name="Schematic description")
    ckt_file_name = models.CharField(max_length=300)
    ckt_sim_case = models.ForeignKey(SimulationCase)


class CircuitSchematicsForm(ModelForm):
    class Meta:
        model = CircuitSchematics
        fields = ('ckt_file_path', 
                'ckt_file_descrip')
