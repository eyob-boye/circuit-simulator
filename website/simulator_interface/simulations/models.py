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


    def clean_sim_time_step(self):
        sim_time_step = float(self.cleaned_data["sim_time_step"])
        if sim_time_step<=0.0:
            raise forms.ValidationError("Simulation time step must be greater than 0.0")
        return sim_time_step

    def clean_sim_time_data(self):
        sim_time_data = float(self.cleaned_data["sim_time_data"])
        if sim_time_data<=0.0:
            raise forms.ValidationError("Data storage rate must be greater than 0.0")
        return sim_time_data

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
        if 'sim_time_step' in cleaned_data and 'sim_time_data' in cleaned_data:
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

    def __unicode__(self):
        return "Circuit spreadsheet "+self.ckt_file_name+" with description "+self.ckt_file_descrip


class CircuitSchematicsForm(ModelForm):
    class Meta:
        model = CircuitSchematics
        fields = ('ckt_file_path', 
                'ckt_file_descrip')


class CircuitComponents(models.Model):
    comp_type = models.CharField(max_length=100)
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50)
    comp_pos = models.CharField(max_length=50)
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200)
    comp_tag = models.CharField(max_length=100)
    sim_case = models.ForeignKey(SimulationCase)

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+self.comp_pos+" in sheet "+self.sheet_name+".csv"


class MeterComponents(models.Model):
    ckt_file_name = models.CharField(max_length=100)
    comp_pos_3D = models.CharField(max_length=50)
    comp_tag = models.CharField(max_length=100)
    comp_type = models.CharField(max_length=25)
    comp_name = models.CharField(max_length=100)
    sim_case = models.ForeignKey(SimulationCase)

    def __unicode__(self):
        return "Meter " + self.comp_type + " with name " + self.comp_tag + \
                " at " + self.comp_pos_3D + " in sheet " + self.ckt_file_name



class ControllableComponents(models.Model):
    ckt_file_name = models.CharField(max_length=100)
    comp_pos_3D = models.CharField(max_length=50)
    comp_tag = models.CharField(max_length=100)
    comp_type = models.CharField(max_length=25)
    comp_name = models.CharField(max_length=100)
    control_tag = models.CharField(max_length=100)
    sim_case = models.ForeignKey(SimulationCase)

    def __unicode__(self):
        return "Comonent " + self.comp_type + " with name " + self.comp_tag + \
                " at " + self.comp_pos_3D + " in sheet " + self.ckt_file_name


class PlotLines(models.Model):
    line_name = models.CharField(max_length=50, verbose_name="Waveform source")
    line_type = models.CharField(max_length=1,
                                        choices=(("M", "Model"),
                                                 ("V", "VariableStorage")),
                                        default="M")
    line_pos = models.IntegerField()
    sim_case = models.ForeignKey(SimulationCase)

    def __unicode__(self):
        return "Source "+self.line_name


class PlotLinesForm(ModelForm):
    class Meta:
        model = PlotLines
        fields = ('line_name', )


class CircuitPlot(models.Model):
    plot_title = models.CharField(max_length=100, default="Plot", \
                verbose_name = "Plot title")
    sim_case = models.ForeignKey(SimulationCase)

    def __unicode__(self):
        return "Plot with title "+self.plot_title


class CircuitPlotForm(ModelForm):
    class Meta:
        model = CircuitPlot
        fields = ('plot_title', )


class CircuitWaveforms(models.Model):
    waveform_legend = models.CharField(max_length=20, blank=True, null=True)
    waveform_scale = models.FloatField(default=1.0, verbose_name = "Scaling Factor")
    circuit_plot = models.ForeignKey(CircuitPlot)
    waveform = models.ManyToManyField(PlotLines)

    def __unicode__(self):
        return "Waveform with legend "+self.waveform_legend


class CircuitWaveformsForm(ModelForm):
    class Meta:
        model = CircuitWaveforms
        fields = ('waveform_legend', \
                'waveform_scale',)


class ControlFile(models.Model):
    control_file = models.FileField(max_length=300, upload_to='', blank=True)
    control_file_name = models.CharField(max_length=300)
    control_file_descrip = models.TextField(blank=True, null=True, verbose_name="Control description")
    sim_case = models.ForeignKey(SimulationCase)

    def __unicode__(self):
        return "Control file with name "+self.control_file_name


class ControlFileForm(ModelForm):
    class Meta:
        model = ControlFile
        fields = ('control_file', \
                'control_file_descrip',)


class ControlInputs(models.Model):
    input_source = models.CharField(max_length=100, \
            verbose_name = "Name of meter")
    input_variable_name = models.CharField(max_length=100, \
            verbose_name = "Desired variable name in control code")
    control_file = models.ForeignKey(ControlFile)

    def __unicode__(self):
        return "Control input with name " + self.input_variable_name + \
                " from source " + self.input_source


class ControlInputsForm(ModelForm):
    class Meta:
        model = ControlInputs
        fields = ('input_variable_name',)


class ControlOutputs(models.Model):
    output_target = models.CharField(max_length=100, \
            verbose_name = "Name of controlled component")
    output_variable_name = models.CharField(max_length=100, \
            verbose_name = "Desired variable name in control code")
    output_initial_value = models.FloatField(verbose_name = "Initial value of output")
    control_file = models.ForeignKey(ControlFile)

    def __unicode__(self):
        return "Control output with name " + self.output_variable_name + \
                " to target " + self.output_target


class ControlOutputsForm(ModelForm):
    class Meta:
        model = ControlOutputs
        fields = ('output_variable_name',\
                'output_initial_value',)


class ControlStaticVariable(models.Model):
    static_variable_name = models.CharField(max_length=100, \
            verbose_name = "Desired variable name in control code")
    static_initial_value = models.FloatField(default=0.0, \
            verbose_name = "Initial value of output")
    control_file = models.ForeignKey(ControlFile)

    def __unicode__(self):
        return "Static variable with name " + self.static_variable_name


class ControlStaticVariableForm(ModelForm):
    class Meta:
        model = ControlStaticVariable
        fields = ('static_variable_name',\
                'static_initial_value',)


class ControlTimeEvent(models.Model):
    time_event_name = models.CharField(max_length=25, \
            verbose_name="Name of time event variable")
    initial_time_value = models.FloatField(verbose_name="Initial time event")
    control_file = models.ForeignKey(ControlFile)

    def __unicode__(self):
        return "Time event variable with name " + self.time_event_name


class ControlTimeEventForm(ModelForm):
    class Meta:
        model = ControlTimeEvent
        fields = ('time_event_name',\
                'initial_time_value',)


class ControlVariableStorage(models.Model):
    variable_storage_name = models.CharField(max_length=100, \
            verbose_name = "Desired variable name in control code")
    storage_initial_value = models.FloatField(verbose_name = "Initial value of variable")
    storage_status = models.CharField(max_length=1,
                                    choices=(("Y", "Yes"),
                                             ("N", "No")),
                                    default="N")
    control_file_name = models.CharField(max_length=100, \
            verbose_name = "Name of control file")
    sim_case = models.ForeignKey(SimulationCase)

    def __unicode__(self):
        return "Variable with name " + self.variable_storage_name


class ControlVariableStorageForm(ModelForm):
    class Meta:
        model = ControlVariableStorage
        fields = ('variable_storage_name',\
                'storage_initial_value',\
                'storage_status',)


class Resistor(models.Model):
    comp_type = models.CharField(max_length=100, default="Resistor", \
            verbose_name="Component type")
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_pos = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200, \
            verbose_name="Found in circuit schematic")
    comp_tag = models.CharField(max_length=100, \
            verbose_name="Component name")
    comp_ckt = models.ForeignKey(CircuitSchematics)

    comp_has_voltage = models.BooleanField(max_length=5, default=False)
    comp_is_meter = models.BooleanField(max_length=5, default=False)
    comp_has_control = models.BooleanField(max_length=5, default=False)
    comp_resistor = models.FloatField(default=100.0, verbose_name="Resistor value")

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+\
                self.comp_pos+" in sheet "+self.sheet_name+\
                ".csv"+" has value "+str(self.comp_resistor)+" Ohms"


class ResistorForm(ModelForm):
    class Meta:
        model = Resistor
        fields = ('comp_resistor', )

    def clean_comp_resistor(self):
        checkres = self.cleaned_data["comp_resistor"]
        if checkres<0.0:
            raise forms.ValidationError("Resistor has to be a positive number.")
        return checkres


class VariableResistor(models.Model):
    comp_type = models.CharField(max_length=100, default="VariableResistor", \
            verbose_name="Component type")
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_pos = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200, \
            verbose_name="Found in circuit schematic")
    comp_tag = models.CharField(max_length=100, \
            verbose_name="Component name")
    comp_ckt = models.ForeignKey(CircuitSchematics)
    
    comp_control_tag = models.CharField(max_length=50, \
            default="Resistance", verbose_name="Control tag")
    comp_control_value = models.FloatField(default=100.0, \
            verbose_name="Controlled resistance")

    comp_has_voltage = models.BooleanField(max_length=5, default=False)
    comp_is_meter = models.BooleanField(max_length=5, default=False)
    comp_has_control = models.BooleanField(max_length=5, default=True)
    comp_resistor = models.FloatField(default=100.0, verbose_name="Initial resistor value")

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+\
                self.comp_pos+" in sheet "+self.sheet_name+\
                ".csv"+" has value "+str(self.comp_resistor)+" Ohms"


class VariableResistorForm(ModelForm):
    class Meta:
        model = VariableResistor
        fields = ('comp_control_tag',
                'comp_resistor',)

    def clean_comp_resistor(self):
        checkres = self.cleaned_data["comp_resistor"]
        if checkres<0.0:
            raise forms.ValidationError("Resistor has to be a positive number.")
        return checkres


class Inductor(models.Model):
    comp_type = models.CharField(max_length=100, default="Inductor", \
            verbose_name="Component type")
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_pos = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200, \
            verbose_name="Found in circuit schematic")
    comp_tag = models.CharField(max_length=100, \
            verbose_name="Component name")
    comp_ckt = models.ForeignKey(CircuitSchematics)

    comp_has_voltage = models.BooleanField(max_length=5, default=False)
    comp_is_meter = models.BooleanField(max_length=5, default=False)
    comp_has_control = models.BooleanField(max_length=5, default=False)
    comp_inductor = models.FloatField(default=0.001, verbose_name="Inductor value")

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+\
                self.comp_pos+" in sheet "+self.sheet_name+\
                ".csv"+" has value "+str(self.comp_inductor)+" Henry"


class InductorForm(ModelForm):
    class Meta:
        model = Inductor
        fields = ('comp_inductor', )

    def clean_comp_inductor(self):
        checkind = self.cleaned_data["comp_inductor"]
        if checkind<0.0:
            raise forms.ValidationError("Inductor has to be a positive number.")
        return checkind


class VariableInductor(models.Model):
    comp_type = models.CharField(max_length=100, default="VariableInductor", \
            verbose_name="Component type")
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_pos = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200, \
            verbose_name="Found in circuit schematic")
    comp_tag = models.CharField(max_length=100, \
            verbose_name="Component name")
    comp_ckt = models.ForeignKey(CircuitSchematics)
    
    comp_control_tag = models.CharField(max_length=50, \
            default="Inductance", verbose_name="Control tag")
    comp_control_value = models.FloatField(default=0.001, \
            verbose_name="Controlled inductance")

    comp_has_voltage = models.BooleanField(max_length=5, default=False)
    comp_is_meter = models.BooleanField(max_length=5, default=False)
    comp_has_control = models.BooleanField(max_length=5, default=True)
    comp_inductor = models.FloatField(default=0.001, verbose_name="Initial inductor value")

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+\
                self.comp_pos+" in sheet "+self.sheet_name+\
                ".csv"+" has value "+str(self.comp_inductor)+" Henry"


class VariableInductorForm(ModelForm):
    class Meta:
        model = VariableInductor
        fields = ('comp_control_tag',
                'comp_inductor',)

    def clean_comp_inductor(self):
        checkind = self.cleaned_data["comp_inductor"]
        if checkind<0.0:
            raise forms.ValidationError("Inductor has to be a positive number.")
        return checkind


class Capacitor(models.Model):
    comp_type = models.CharField(max_length=100, default="Capacitor", \
            verbose_name="Component type")
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_pos = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200, \
            verbose_name="Found in circuit schematic")
    comp_tag = models.CharField(max_length=100, \
            verbose_name="Component name")
    comp_ckt = models.ForeignKey(CircuitSchematics)
    
    comp_polarity_3D = models.CharField(max_length=50, \
            verbose_name="Positive polarity towards")
    comp_polarity = models.CharField(max_length=50, \
            verbose_name="Positive polarity towards")

    comp_has_voltage = models.BooleanField(max_length=5, default=True)
    comp_is_meter = models.BooleanField(max_length=5, default=False)
    comp_has_control = models.BooleanField(max_length=5, default=False)
    comp_capacitor = models.FloatField(default=10.0e-6, verbose_name="Capacitor value")

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+\
                self.comp_pos+" in sheet "+self.sheet_name+".csv"+\
                " has value "+str(self.comp_capacitor)+" Farad"


class CapacitorForm(ModelForm):
    class Meta:
        model = Capacitor
        fields = ('comp_capacitor', \
                'comp_polarity')

    def clean_comp_capacitor(self):
        checkcap = self.cleaned_data["comp_capacitor"]
        if checkcap<0.0:
            raise forms.ValidationError("Capcitor has to be a positive number.")
        return checkcap


class Voltage_Source(models.Model):
    comp_type = models.CharField(max_length=100, default="VoltageSource", \
            verbose_name="Component type")
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_pos = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200, \
            verbose_name="Found in circuit schematic")
    comp_tag = models.CharField(max_length=100, \
            verbose_name="Component name")
    comp_ckt = models.ForeignKey(CircuitSchematics)
    
    comp_polarity_3D = models.CharField(max_length=50, \
            verbose_name="Positive polarity towards")
    comp_polarity = models.CharField(max_length=50, \
            verbose_name="Positive polarity towards")

    comp_has_voltage = models.BooleanField(max_length=5, default=True)
    comp_is_meter = models.BooleanField(max_length=5, default=False)
    comp_has_control = models.BooleanField(max_length=5, default=False)
    comp_volt_peak = models.FloatField(default=120.0, verbose_name="Peak voltage")
    comp_volt_freq = models.FloatField(default=60.0, verbose_name="Voltage frequency")
    comp_volt_phase = models.FloatField(default=0.0, verbose_name="Phase angle (degrees)")
    comp_volt_offset = models.FloatField(default=0.0, verbose_name="Dc offset")

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+\
                self.comp_pos+" in sheet "+self.sheet_name+".csv"+\
                " has peak value "+str(self.comp_volt_peak)+" Volts"


class Voltage_SourceForm(ModelForm):
    class Meta:
        model = Voltage_Source
        fields = ('comp_volt_peak', \
                'comp_volt_peak', \
                'comp_volt_freq', \
                'comp_volt_phase', \
                'comp_volt_offset', \
                'comp_polarity')

    def clean_comp_volt_peak(self):
        checkcomp_volt_peak = self.cleaned_data["comp_volt_peak"]
        if checkcomp_volt_peak<0.0:
            raise forms.ValidationError("Peak voltage has to be a positive number.")
        return checkcomp_volt_peak

    def clean_comp_volt_freq(self):
        checkcomp_volt_freq = self.cleaned_data["comp_volt_freq"]
        if checkcomp_volt_freq<0.0:
            raise forms.ValidationError("Frequency has to be a positive number.")
        return checkcomp_volt_freq


class Controlled_Voltage_Source(models.Model):
    comp_type = models.CharField(max_length=100, default="ControlledVoltageSource", \
            verbose_name="Component type")
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_pos = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200, \
            verbose_name="Found in circuit schematic")
    comp_tag = models.CharField(max_length=100, \
            verbose_name="Component name")
    comp_ckt = models.ForeignKey(CircuitSchematics)
    
    comp_control_tag = models.CharField(max_length=50, \
            default="Voltage", verbose_name="Control tag")
    comp_control_value = models.FloatField(default=0.0, \
            verbose_name="Controlled Voltage")
    
    comp_polarity_3D = models.CharField(max_length=50, \
            verbose_name="Positive polarity towards")
    comp_polarity = models.CharField(max_length=50, \
            verbose_name="Positive polarity towards")

    comp_has_voltage = models.BooleanField(max_length=5, default=True)
    comp_is_meter = models.BooleanField(max_length=5, default=False)
    comp_has_control = models.BooleanField(max_length=5, default=True)
    comp_voltage = models.FloatField(default=0.0, verbose_name="Initial voltage")

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+\
                self.comp_pos+" in sheet "+self.sheet_name+".csv"+\
                " has peak value "+str(self.comp_volt_peak)+" Volts"


class Controlled_Voltage_SourceForm(ModelForm):
    class Meta:
        model = Controlled_Voltage_Source
        fields = ('comp_control_tag', \
                'comp_control_value',)


class Ammeter(models.Model):
    comp_type = models.CharField(max_length=100, default="Ammeter", \
            verbose_name="Component type")
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_pos = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200, \
            verbose_name="Found in circuit schematic")
    comp_tag = models.CharField(max_length=100, \
            verbose_name="Component name")
    comp_ckt = models.ForeignKey(CircuitSchematics)
    
    comp_polarity_3D = models.CharField(max_length=50, \
            verbose_name="Positive polarity towards")
    comp_polarity = models.CharField(max_length=50, \
            verbose_name="Positive polarity towards")

    comp_has_voltage = models.BooleanField(max_length=5, default=False)
    comp_is_meter = models.BooleanField(max_length=5, default=True)
    comp_has_control = models.BooleanField(max_length=5, default=False)

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+\
                self.comp_pos+" in sheet "+self.sheet_name+".csv"


class AmmeterForm(ModelForm):
    class Meta:
        model = Ammeter
        fields = ('comp_polarity', )


class Voltmeter(models.Model):
    comp_type = models.CharField(max_length=100, default="Voltmeter", \
            verbose_name="Component type")
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_pos = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200, \
            verbose_name="Found in circuit schematic")
    comp_tag = models.CharField(max_length=100, \
            verbose_name="Component name")
    comp_ckt = models.ForeignKey(CircuitSchematics)
    
    comp_polarity_3D = models.CharField(max_length=50, \
            verbose_name="Positive polarity towards")
    comp_polarity = models.CharField(max_length=50, \
            verbose_name="Positive polarity towards")
    
    comp_volt_level = models.FloatField(default=120.0, \
            verbose_name="Rated voltage level to be measured")

    comp_has_voltage = models.BooleanField(max_length=5, default=False)
    comp_is_meter = models.BooleanField(max_length=5, default=True)
    comp_has_control = models.BooleanField(max_length=5, default=False)

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+\
                self.comp_pos+" in sheet "+self.sheet_name+".csv"


class VoltmeterForm(ModelForm):
    class Meta:
        model = Voltmeter
        fields = ('comp_volt_level', \
                'comp_polarity')

    def clean_comp_volt_level(self):
        checkcomp_volt_level = self.cleaned_data["comp_volt_level"]
        if checkcomp_volt_level<0.0:
            raise forms.ValidationError("Voltage level has to be a positive number.")
        return checkcomp_volt_level


class Diode(models.Model):
    comp_type = models.CharField(max_length=100, default="Diode", \
            verbose_name="Component type")
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_pos = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200, \
            verbose_name="Found in circuit schematic")
    comp_tag = models.CharField(max_length=100, \
            verbose_name="Component name")
    comp_ckt = models.ForeignKey(CircuitSchematics)
    
    comp_polarity_3D = models.CharField(max_length=50, \
            verbose_name="Cathode polarity towards")
    comp_polarity = models.CharField(max_length=50, \
            verbose_name="Cathode polarity towards")
    
    comp_volt_level = models.FloatField(default=120.0, \
            verbose_name="Rated voltage level")

    comp_has_voltage = models.BooleanField(max_length=5, default=True)
    comp_is_meter = models.BooleanField(max_length=5, default=False)
    comp_has_control = models.BooleanField(max_length=5, default=False)

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+\
                self.comp_pos+" in sheet "+self.sheet_name+".csv"


class DiodeForm(ModelForm):
    class Meta:
        model = Diode
        fields = ('comp_volt_level', \
                'comp_polarity')

    def clean_comp_volt_level(self):
        checkcomp_volt_level = self.cleaned_data["comp_volt_level"]
        if checkcomp_volt_level<0.0:
            raise forms.ValidationError("Voltage level has to be a positive number.")
        return checkcomp_volt_level


class Switch(models.Model):
    comp_type = models.CharField(max_length=100, default="Switch", \
            verbose_name="Component type")
    comp_number = models.IntegerField()
    comp_pos_3D = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_pos = models.CharField(max_length=50, \
            verbose_name="Component position")
    comp_sheet = models.IntegerField()
    sheet_name = models.CharField(max_length=200, \
            verbose_name="Found in circuit schematic")
    comp_tag = models.CharField(max_length=100, \
            verbose_name="Component name")
    comp_ckt = models.ForeignKey(CircuitSchematics)
    
    comp_polarity_3D = models.CharField(max_length=50, \
            verbose_name="Negative polarity towards")
    comp_polarity = models.CharField(max_length=50, \
            verbose_name="Negative polarity towards")

    comp_control_tag = models.CharField(max_length=50, \
            default="Gate", verbose_name="Control tag")
    comp_control_value = models.FloatField(default=0.0, \
            verbose_name="Gate signal")
    
    comp_volt_level = models.FloatField(default=120.0, \
            verbose_name="Rated voltage level")

    comp_has_voltage = models.BooleanField(max_length=5, default=True)
    comp_is_meter = models.BooleanField(max_length=5, default=False)
    comp_has_control = models.BooleanField(max_length=5, default=True)

    def __unicode__(self):
        return "Component "+self.comp_type+" with name "+self.comp_tag+" at "+\
                self.comp_pos+" in sheet "+self.sheet_name+".csv"


class SwitchForm(ModelForm):
    class Meta:
        model = Switch
        fields = ('comp_volt_level', \
                'comp_polarity', \
                'comp_control_tag',)

    def clean_comp_volt_level(self):
        checkcomp_volt_level = self.cleaned_data["comp_volt_level"]
        if checkcomp_volt_level<0.0:
            raise forms.ValidationError("Voltage level has to be a positive number.")
        return checkcomp_volt_level

