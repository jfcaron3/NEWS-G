#!/usr/bin/python
# This is an adaptation of the original DD2Estimator Jupyter/IPython
# notebook by Paco (Francisco Vazquez de Sola Fernandez <14favd@queensu.ca>)
# Adaptation by Jean-Francois Caron (jfc3@queensu.ca) from a copy made in 
# January 2021. 
# This version is a standalone python program and can be run directly.
# You should edit the path for the libquadis.so library, the default 
# input files, and the output plot folder.
# You can also import it as a module to use the DDFilter function yourself.
# Note that this adaptation ignores the specialized treatment of simulated
# data files.
#
# Usage:
# python PulseExplorer.py [T1Filename.root] [T2Filename.root] [jsonfilename]
#
# This will draw pulses from the T1 file, processed using parameters from 
# the JSON file.  Cuts are applied using variables from the T2 file.
#
# NOTE: Make sure the JSON file you enter is the same as was used to 
#       generate the T2 file, otherwise this won't make any sense!
# Note: If you enter an invalid T2Filename (e.g. "None"), then T2 cuts are
#       not available. This way you can use Pulse Explorer without a T2 file.
#       You can still enter cuts on the T1 variable names, but these are
#       very limited since most data is inside the NEvent.

from __future__ import print_function, division
import ROOT
import numpy
import os, time, sys

# For use on JF's system while developping
ROOT.gSystem.Load("/home/jfcaron/Quadis/quadis/build/lib/libquadis.so")
default_datadir = "/home/jfcaron/NEWS-G/RMTL Test 6/SPC Data/"
default_T1 = os.path.join(default_datadir,"ul03x001_T1.root")
default_T2 = os.path.join(default_datadir,"ul03x001_q00.root")
default_json = os.path.join(default_datadir,"ul03x001_q00_input.json")

# For use on sphere server
#ROOT.gSystem.Load('/usr/local/software/quadis/latest/quadis/build/lib/libquadis.so')
#default_datadir = "/home/newsdata/scripts"
#default_T1 = os.path.join(default_datadir,"PulseExplorer_example_T1.root")
#default_T2 = os.path.join(default_datadir,"PulseExplorer_example_T2.root")
#default_json = os.path.join(default_datadir,"PulseExplorer_example.json")

output_dir = os.path.join(default_datadir,"plots")
output_formats = ['png'] # Add 'pdf', 'root' etc for other/more formats

# Only use the defaults if the user did not provide input files.
try:
    T1_filename = sys.argv[1]
except IndexError:
    T1_filename = default_T1
try:
    T2_filename = sys.argv[2]
except IndexError:
    T2_filename = default_T2
try:
    json_filename = sys.argv[3]
except IndexError:
    json_filename = default_json

runname = os.path.splitext(os.path.basename(T1_filename))[0]

# Open the T1 file.
T1_file = ROOT.TFile(T1_filename,"READ")
T1_tree = T1_file.Get("T1")
nEntries = T1_tree.GetEntries()
print("The T1 tree has {:,} entries.".format(nEntries))

# Open the T2 file.
T2_file = ROOT.TFile(T2_filename,"READ")
if T2_file.IsZombie():
    print("T2 is not a valid file, T2 cuts not available.")
    T2 = False
else:
    T2 = True
    T2_tree = T2_file.Get("T2")
    T2_nEntries = T2_tree.GetEntries()
    print("The T2 tree has {:,} entries.".format(T2_nEntries))

# Create an NEvent object to store the retrieved events.
nevent = ROOT.NEvent()
T1_tree.SetBranchAddress("NEvent",nevent)

# Print the sampling frequency, assuming it does not change during 
# the run or between channels.
T1_tree.GetEntry(0)
npulse = nevent.GetPulse(0)
SamplingPeriod = nevent.GetPulse(0).GetSamplingPeriod()
print("Sampling frequency is {:,.6f} MHz\n".format(1000./SamplingPeriod))

# Calculate the duration of the run and the average event rate.
T1_tree.GetEntry(0) # Get the first event
GlobalTime0 = nevent.GetPulseData(0).GetParameterByName("GlobalTime")
T1_tree.GetEntry(nEntries-1) # Get the last event
GlobalTimeN = nevent.GetPulseData(0).GetParameterByName("GlobalTime")
runhours = (GlobalTimeN-GlobalTime0)/3600.0
print("Total length of the run is {:,.2f} h".format(runhours))
print("The event rate is {:,.2f} Hz".format(nEntries/runhours/3600.0) )
print("\n---------------")

# Print out all parameter names for the pulse data
num_pulses = nevent.GetNPulses()
n = nevent.GetPulseData(0).GetNParameters()
print("The pulse data has "+str(n)+" parameters of type double available")
paramList = nevent.GetPulseData(0).GetParameterList()
for i in range(n):
    print(str(i)+": "+paramList.GetParameterName(i))
print()
n = nevent.GetPulseData(0).GetNParametersUL()
print("The pulse data has "+str(n)+" parameters of type size_t available")
paramList = nevent.GetPulseData(0).GetParameterListUL()
for i in range(n):
    print(str(i)+": "+paramList.GetParameterName(i))

# This is a mega-object from Quadis that does most of the processing.
ddFilter = ROOT.NFDDMethod() 
ddFilter.SetFFToption(0) #0-3, trade start up time for processing time 

# Create an NPulse object to hold retrieved data.
ddPulse = ROOT.NPulse() # Will hold the double-deconvolved pulse.

class Quiet:
    """Context manager for silencing certain ROOT operations.  Usage:
    with Quiet(level = ROOT.kInfo+1):
       foo_that_makes_output

    You can set a higher or lower warning level to ignore different
    kinds of messages.  After the end of indentation, the level is set
    back to what it was previously.
    """
    def __init__(self, level=ROOT.kInfo + 1):
        self.level = level

    def __enter__(self):
        self.oldlevel = ROOT.gErrorIgnoreLevel
        ROOT.gErrorIgnoreLevel = self.level

    def __exit__(self, type, value, traceback):
        ROOT.gErrorIgnoreLevel = self.oldlevel

def DDFilter(OutputPulseType, EventNumber, PulseNumber, InLogPath, xunit):
    # Get the requested Entry, Pulse and the PulseData
    T1_tree.GetEntry(EventNumber)
    npulse = nevent.GetPulse(int(PulseNumber))
    ndata = nevent.GetPulseData(int(PulseNumber))
    sP = npulse.GetSamplingPeriod() # In milliseconds.

    # If the user requested converting to microseconds or stay in samples.
    conversionfactor = 1000./sP # This is the sampling frequency in MHz.
    if xunit == "microseconds":
        samplingconversion = 1. # Keep output units in microseconds
    else:
        samplingconversion = conversionfactor # Convert output units to samples.

    SC = samplingconversion
    CF = conversionfactor

    # Something to do with a trapezoidal filter
    # fTrapGapSamples + fTrapRiseSamples/2 + 1 from NFDDMethod.hh in Quadis
    trapcorr = ddFilter.GetTrapGapSamples() + ddFilter.GetTrapRiseSamples()//2+1

    # Set the type of output pulse you want to plot.
    ddFilter.SetOutputPulseType(OutputPulseType)

    # I ignored the part of the code where you can set individual ddFilter
    # variables in python, it only reads from the 'JSON' file.
    ddFilter.ReadJSONInput(InLogPath)

    # Run the mega-algorithm.  The output is put in ddPulse
    ddFilter.Filter(npulse,ddPulse)

    # Get output parameters from the mega-algorithm into python variables. 
    # I hate that this changes some of the names from Quadis.
    DDTriggerTime = ddFilter.GetCalcParameterByName("StartTime")
    DDBaseline = ddFilter.GetCalcParameterByName("Baseline")
    DDSmoothBaseline = ddFilter.GetCalcParameterByName("Smooth5Baseline")
    DDSmoothBaseline_RMS = ddFilter.GetCalcParameterByName("Smooth5Baseline_RMS")
    DDAmplitude = ddFilter.GetCalcParameterByName("RawAmpl")
    DDMaxSample = ddFilter.GetCalcParameterByName("RawMaxSample")
    DD10pctFMR = ddFilter.GetCalcParameterByName("RawRise10pct")
    DD90pctFMR = ddFilter.GetCalcParameterByName("RawRise90pct")
    DDRawRise = ddFilter.GetCalcParameterByName("RawRise")
    DD50pctFMR = ddFilter.GetCalcParameterByName("RawRise50pct")
    DD50pctFMF = ddFilter.GetCalcParameterByName("RawFall50pct")
    DDWidth = ddFilter.GetCalcParameterByName("RawWidth")
    DDStopTime = ddFilter.GetCalcParameterByName("StopTime")
    DDThresholdStop = ddFilter.GetCalcParameterByName("ThresholdStop")
    DDDecMinimum = ddFilter.GetCalcParameterByName("DecMinimum")
    DDDecBaselineSlopeStart = ddFilter.GetCalcParameterByName("DecBaselineSlope_start")
    DDDecBaselineSlopeEnd = ddFilter.GetCalcParameterByName("DecBaselineSlope_end")
    DDRealTriggerTime = ddFilter.GetCalcParameterByName("RealStartTime")
    DDRealStopTime = ddFilter.GetCalcParameterByName("RealStopTime")
    DDDecBaselineStart = ddFilter.GetCalcParameterByName("BaselineStart")
    DDDecBaselineEnd = ddFilter.GetCalcParameterByName("BaselineEnd")
    DDAmpl = ddFilter.GetCalcParameterByName("Ampl")
    DDAmplADU = ddFilter.GetCalcParameterByName("AmplADU")
    DDRise = ddFilter.GetCalcParameterByName("Rise")
    DD10pct = ddFilter.GetCalcParameterByName("Rise10pct")
    DD25pct = ddFilter.GetCalcParameterByName("Rise25pct")
    DD75pct = ddFilter.GetCalcParameterByName("Rise75pct")
    DD90pct = ddFilter.GetCalcParameterByName("Rise90pct")

    # I ignored the parts for variables of simulated data and peak-finding.

    SambaAmpl = ndata.GetParameterByName("Amplitude")
    SambaRise = ndata.GetParameterByName("Risetime")
    SambaBaseline = ndata.GetParameterByName("Baseline")
    SambaBaselineNoise = ndata.GetParameterByName("BaselineNoise")
    DeltaTPrevS = ndata.GetParameterByName("DeltaTPrevS")

    # More variables that were passed into the function in DD2Estimator
    # as arguments and SET into the ddFilter, now they are already there
    # from having read the JSON file, so we have to GET them.
    SamplesFromBeforeTrigger = ddFilter.GetSamplesFromBeforeTrigger()
    RawBaselineSamples = ddFilter.GetRawBaselineSamples()
    StartSafetyMargin  = ddFilter.GetStartSafetyMargin()
    StartThreshold = ddFilter.GetThresholdStart() # Note: words reversed.
    StopThreshold = ddFilter.GetThresholdStop() # Note: words reversed.
    DecBaselineSamples = ddFilter.GetDecBaselineSamples()
    DecBaselineSlopeSamples = ddFilter.GetDecBaselineSlopeSamples()
    EndSafetyMargin = ddFilter.GetEndSafetyMargin()
    Charge2Volt = ddFilter.GetCharge2Volt()
    Volt2ADU = ddFilter.GetVolt2ADU()

    if DDTriggerTime<0 and DDMaxSample<0 and DDWidth<0:
        print("Warning: Empty/unphysical pulse.\n")

    # Here DD2Estimator would print all the variables to the terminal.
    # You can print them yourself if you want.

    # Names and Y-axis labels for the 11 different pulse types.
    pulsenames = {-1:"Raw Pulse",0:"Baseline Removed",1:"Trapezoidal Filter",
                   2:"Baseline Removed, Light Smoothing",
                   3:"Deconvolved from Preamp",
                   4:"Deconvolved from Preamp, Smoothed",
                   5:"Double Deconvolved, Harsh Smoothing",
                   6:"Double Deconvolved, Light Smoothing",
                   7:"Integral of Double Deconvolved Pulse",
                   8:"Baseline Removed Integral of Double Deconvolved Pulse",
                   9:"Baseline Removed Integral of Double Deconvolved Pulse (in ADUs)"}
    yaxislabels = {-1:"ADUs",0:"ADUs",1:"ADUs",2:"ADUs",
                  3:"Electrons per ns",4:"Electrons per ns",
                  5:"Electrons per ns",6:"Electrons per ns",
                  7:"Electrons",8:"Electrons",9:"ADUs"}

    # Prepare some arrays with the pulse and processed pulse data.
    if OutputPulseType==-1: # "raw pulse"
        # Get the values from the INPUT pulse.
        dd_array = numpy.array( npulse.GetArray() )
    elif OutputPulseType==1:# " trapezoidal filter"
        # Get the 1th element from the OUTPUT pulse.
        trapcorr0 = ddPulse.GetArray()[1]
        # Turn the OUTPUT pulse into a list?
        preddarray = list(ddPulse.GetArray())
        # Recall trapcorr is the combo of trapezoidal filter parameters.
        # So this is some kind of re-arranging or scaling of the output.
        dd_array = numpy.array( [trapcorr0]*trapcorr + preddarray[:-trapcorr] )
    else: # All non-raw and non-trapezoidal filter output types.
        dd_array = numpy.array( ddPulse.GetArray() )
    
    if OutputPulseType == 9: #Change unit from electrons to ADUs
        dd_array = dd_array * (Charge2Volt*Volt2ADU)
    
    # Make a TGraph to hold the pulse data, the actual X vs Y waveform.
    N = len(dd_array)
    endpoint = N*SC/CF
    g_x = numpy.linspace(0, endpoint, N, endpoint = False, dtype=numpy.double)
    g_y = numpy.array(dd_array,dtype=numpy.double)
    g = ROOT.TGraph(N,g_x,g_y)
    g.SetLineColor(ROOT.kAzure)
    g.SetMarkerColor(ROOT.kAzure)
    g.Draw("ALP")
    g.SetTitle("%s, %d, %d" % (runname,EventNumber,PulseNumber))
    g.SetName(pulsenames[OutputPulseType])
    g.GetYaxis().SetTitle(yaxislabels[OutputPulseType])
    g.GetXaxis().SetTitle(xunit)
    
    # Now I make all the markup lines for the different types of plots.
    narrow = 3
    wide = 7
    tl_red = ROOT.TLine()
    tl_red.SetLineWidth(narrow)
    tl_red.SetLineColor(ROOT.kRed)

    tl_red_wide = ROOT.TLine()
    tl_red_wide.SetLineWidth(wide)
    tl_red_wide.SetLineColor(ROOT.kRed)

    tl_green = ROOT.TLine()
    tl_green.SetLineWidth(narrow)
    tl_green.SetLineColor(ROOT.kGreen)

    tl_green_wide = ROOT.TLine()
    tl_green_wide.SetLineWidth(wide)
    tl_green_wide.SetLineColor(ROOT.kGreen)

    tl_green_dotted = ROOT.TLine()
    tl_green_dotted.SetLineWidth(narrow)
    tl_green_dotted.SetLineColor(ROOT.kGreen+2)
    tl_green_dotted.SetLineStyle(ROOT.kDotted)

    tl_green_dashed = ROOT.TLine()
    tl_green_dashed.SetLineWidth(narrow)
    tl_green_dashed.SetLineColor(ROOT.kGreen+2)
    tl_green_dashed.SetLineStyle(ROOT.kDashed)

    tl_green_dotdashed = ROOT.TLine()
    tl_green_dotdashed.SetLineWidth(narrow)
    tl_green_dotdashed.SetLineColor(ROOT.kGreen+2)
    tl_green_dotdashed.SetLineStyle(ROOT.kDashDotted)

    tl_black = ROOT.TLine()
    tl_black.SetLineWidth(narrow)
    tl_black.SetLineColor(ROOT.kBlack)

    tl_black_dashed = ROOT.TLine()
    tl_black_dashed.SetLineWidth(narrow)
    tl_black_dashed.SetLineColor(ROOT.kBlack)
    tl_black_dashed.SetLineStyle(ROOT.kDashed)

    tl_black_dotdashed = ROOT.TLine()
    tl_black_dotdashed.SetLineWidth(narrow)
    tl_black_dotdashed.SetLineColor(ROOT.kBlack)
    tl_black_dotdashed.SetLineStyle(ROOT.kDashDotted)

    tl_darkgray = ROOT.TLine()
    tl_darkgray.SetLineWidth(narrow)
    tl_darkgray.SetLineColor(ROOT.kGray)

    tl_darkgray_dashed = ROOT.TLine()
    tl_darkgray_dashed.SetLineWidth(narrow)
    tl_darkgray_dashed.SetLineColor(ROOT.kGray)
    tl_darkgray_dashed.SetLineStyle(ROOT.kDashed)

    te_red = ROOT.TEllipse()
    te_red.SetLineColor(ROOT.kRed)
    te_red.SetLineWidth(narrow)
    te_red.SetFillStyle(0)
    
    te_green = ROOT.TEllipse()
    te_green.SetLineColor(ROOT.kGreen)
    te_green.SetLineWidth(narrow)
    te_green.SetFillStyle(0)
    
    legend = ROOT.TLegend(0.1,0.9,0.3,0.6)
    legend.AddEntry(g,pulsenames[OutputPulseType])

    if OutputPulseType in [-1]: #Raw pulse
        # "raw1" Window over which baseline is calculated.
        # DDTriggerTime is StartTime in Quadis
        # SamplesFromBeforeTrigger: number of samples before the trigger 
        #                           for calculating baseline
        # RawBaselineSamples: number of samples to use to calculate baseline
        # So x2_raw1 is the time where the baseline calculation starts.
        #    x1_raw1 is the time where the baseline calculation ends.
        x1_raw1 = (DDTriggerTime-SamplesFromBeforeTrigger-RawBaselineSamples)*SC/CF
        x2_raw1 = (DDTriggerTime-SamplesFromBeforeTrigger)*SC/CF
        y1_raw1 = y2_raw1 = DDBaseline
        tl_red.DrawLine(x1_raw1,y1_raw1,x2_raw1,y2_raw1)
        legend.AddEntry(tl_red,"Baseline Window","l")
        # "raw2"
        # In the original DD2Estimator "raw2" was a bigger red line to help
        # the user find the baseline window, cuz "raw1" is typically pretty 
        # short in the initial view. It confusingly used the 
        # StartSafetyMargin variable to set where to draw.
        # I changed it to an ellipse around the "raw1" line.
        xax = g.GetXaxis()
        radius = abs(xax.GetXmax() - xax.GetXmin())/30.0
        te_red.DrawEllipse(x1_raw1, y1_raw1, radius, radius, 0, 360,0)
    if OutputPulseType in [0,2]: #Baseline removed / smoothed
        # "br1" Vertical line at the 50% point on the rising edge.
        x1_br1 = x2_br1 = DD50pctFMR*SC # n.b.: This is RawRise50pct
        y1_br1, y2_br1 = [0,DDAmplitude]
        tl_green.DrawLine(x1_br1,y1_br1,x2_br1,y2_br1)
        legend.AddEntry(tl_green,"50% Rising Amplitude","l")
        # "br2" Vertical line at the 50% point on the falling edge.
        x1_br2 = x2_br2 = DD50pctFMF*SC # n.b.: This is RawFall50pct
        y1_br2, y2_br2 = [0,DDAmplitude]
        tl_green.DrawLine(x1_br2,y1_br2,x2_br2,y2_br2)
        legend.AddEntry(tl_green,"50% Falling Amplitude","l")
        # "br3" Horizontal line joining br1 and br2 midpoints.
        x1_br3, x2_br3 = [DD50pctFMR*SC, DD50pctFMF*SC]
        y1_br3 = y2_br3 = 0.5*DDAmplitude
        tl_green_dotted.DrawLine(x1_br3,y1_br3,x2_br3,y2_br3)
        legend.AddEntry(tl_green_dotted,"Connecting 50% Points","l")
        # "br4" Vertical line indicating maximum point of pulse.
        x1_br4 = x2_br4 = DDMaxSample*SC/CF # n.b.: This is RawMaxSample
        y1_br4, y2_br4 = [0,1.05*DDAmplitude]
        tl_green_dashed.DrawLine(x1_br4,y1_br4,x2_br4,y2_br4)
        legend.AddEntry(tl_green_dashed,"Maximum Amplitude","l")
        # "br5" Vertical line indicating the 10% point on the rising edge
        x1_br5 = x2_br5 = DD10pctFMR*SC # 10% rising edge
        y1_br5, y2_br5 = [0,DDAmplitude]
        tl_black_dotdashed.DrawLine(x1_br5,y1_br5,x2_br5,y2_br5)
        legend.AddEntry(tl_green_dotdashed,"10% Rising Amplitude","l")
        # "br6" Vertical line indicating the 90% point on the rising edge
        x1_br6 = x2_br6 = DD90pctFMR*SC
        y1_br6, y2_br6 = [0,DDAmplitude]
        tl_black_dotdashed.DrawLine(x1_br6,y1_br6,x2_br6,y2_br6)
        legend.AddEntry(tl_black_dotdashed,"90% Rising Amplitude","l")
        # "br7" Window over which baseline is calculated.
        # This is the same as "raw1" but with y = 0 instead of y = Baseline
        x1_br7 = (DDTriggerTime-SamplesFromBeforeTrigger-RawBaselineSamples)*SC/CF
        x2_br7 = (DDTriggerTime-SamplesFromBeforeTrigger)*SC/CF
        y1_br7 = y2_br7 = 0
        tl_red_wide.DrawLine(x1_br7,y1_br7,x2_br7,y2_br7)
        legend.AddEntry(tl_red_wide,"Baseline Window","l")
        # "br8" is like "raw2" from the -1 case above.
        # I draw an ellipse instead.
        xax = g.GetXaxis()
        radius = abs(xax.GetXmax() - xax.GetXmin())/30.0
        te_red.DrawEllipse(x1_br7, y1_br7, radius, radius, 0, 360,0)
    elif OutputPulseType in [1]: #Trapezoidal filter
        # "tr1" Vertical line where the trigger level was crossed.
        x1_tr1 = x2_tr1 = DDTriggerTime*SC/CF
        y1_tr1, y2_tr1 = [0,DDAmplitude]
        tl_black_dashed.DrawLine(x1_tr1,y1_tr1,x2_tr1,y2_tr1)
        legend.AddEntry(tl_black_dashed,"Trigger Crossing Time","l")
        # tr2source
        x1_tr2 = DDTriggerTime*SC/CF
        x2_tr2 = DDStopTime*SC/CF # n.b.: This is StopTime
        y1_tr2 = y2_tr2 = StartThreshold
        tl_red.DrawLine(x1_tr2,y1_tr2,x2_tr2,y2_tr2)
        legend.AddEntry(tl_red,"Trigger Crossing Amplitude","l")
    elif OutputPulseType in [5]: #Double deconvolved, harsh smoothing
        # "de1" Vertical line at DDRealTriggerTime
        # DDRealTriggerTime is the trigger time minus a fixed SafetyMargin.
        x1_de1 = x2_de1 = DDRealTriggerTime*SC # n.b.: This is RealStartTime
        y1_de1, y2_de1 = [0,DDAmplitude/25]
        tl_darkgray_dashed.DrawLine(x1_de1,y1_de1,x2_de1,y2_de1)
        legend.AddEntry(tl_darkgray_dashed,"DDRealTriggerTime (Integration Start)","l")
        # "de2" Vertical line at DDStopTime
        # DDStopTime is the time at which the pulse drops below the
        # the threshold.
        x1_de2 = x2_de2 = DDStopTime*SC/CF
        y1_de2, y2_de2 = [0,DDAmplitude/25]
        tl_black_dashed.DrawLine(x1_de2,y1_de2,x2_de2,y2_de2)
        legend.AddEntry(tl_black_dashed,"DDStopTime (Integration Stop)","l")
        # "de3" Horizontal line at StopThreshold
        # StopThreshold is is the threshold argument given to DDFilter
        x1_de3, x2_de3 = [DDRealTriggerTime*SC,DDStopTime*SC/CF]
        y1_de3 = y2_de3 = StopThreshold
        tl_darkgray.DrawLine(x1_de3,y1_de3,x2_de3,y2_de3)
        legend.AddEntry(tl_darkgray,"StopThreshold (Threshold for Integration Stop)","l")
        # "de4" Horizontal line at DDThresholdStop
        # DDThresholdStop is the actual value used if ThresholdStop/StopThreshold
        # wasn't crossed.
        x1_de4, x2_de4 = [DDRealTriggerTime*SC,DDStopTime*SC/CF]
        y1_de4 = y2_de4 = DDThresholdStop
        tl_black.DrawLine(x1_de4,y1_de4,x2_de4,y2_de4)
        legend.AddEntry(tl_black,"DDThresholdStop (Value Used for Integration Stop)","l")
    elif OutputPulseType in [7]: #Integral of deconvolved
        # These shifted Y-coordinates are used several times here.
        # They are just convenient places to start/end vertical lines.
        Y1_opt7 = -0.05*DDAmpl+DDDecBaselineStart
        Y2_opt7 =  1.05*DDAmpl+DDDecBaselineStart
        # id1source
        x1_id1 = x2_id1 = DD10pct*SC
        y1_id1, y2_id1 = Y1_opt7, Y2_opt7
        tl_green.DrawLine(x1_id1,y1_id1,x2_id1,y2_id1)
        legend.AddEntry(tl_green,"10% Amplitude","l")
        #id1bsource
        x1_id1b = x2_id1b = DD25pct*SC
        y1_id1b, y2_id1b = Y1_opt7, Y2_opt7
        tl_green_dashed.DrawLine(x1_id1b,y1_id1b,x2_id1b,y2_id1b)
        legend.AddEntry(tl_green_dashed,"25% Amplitude","l")
        #id2bsource
        x1_id2b = x2_id2b = DD75pct*SC
        y1_id2b, y2_id2b = Y1_opt7, Y2_opt7
        tl_green_dashed.DrawLine(x1_id2b,y1_id2b,x2_id2b,y2_id2b)
        legend.AddEntry(tl_green_dashed,"75% Amplitude","l")
        # id2source
        x1_id2 = x2_id2 = DD90pct*SC
        y1_id2, y2_id2 = Y1_opt7, Y2_opt7
        tl_green.DrawLine(x1_id2,y1_id2,x2_id2,y2_id2)
        legend.AddEntry(tl_green,"90% Amplitude","l")
        #id3source
        x1_id3 = x2_id3 = DDRealStopTime*SC
        y1_id3, y2_id3 = Y1_opt7, Y2_opt7
        tl_darkgray_dashed.DrawLine(x1_id3,y1_id3,x2_id3,y2_id3)
        legend.AddEntry(tl_darkgray_dashed,"Integration Stop Time","l")
        #id4source
        x1_id4 = x2_id4 = DDRealTriggerTime*SC
        y1_id4, y2_id4 = Y1_opt7, Y2_opt7
        tl_darkgray_dashed.DrawLine(x1_id4,y1_id4,x2_id4,y2_id4)
        legend.AddEntry(tl_darkgray_dashed,"Integration Start Time","l")
        #id7source
        x1_id7 = (DDRealTriggerTime - DecBaselineSamples/CF)*SC
        x2_id7 = DDRealTriggerTime*SC
        y1_id7 = y2_id7 = DDDecBaselineStart
        tl_green_wide.DrawLine(x1_id7,y1_id7,x2_id7,y2_id7)
        legend.AddEntry(tl_green_wide,"Baseline Before Pulse","l")
        # id8source
        x1_id8 = DDRealStopTime*SC
        x2_id8 = (DDRealStopTime + DecBaselineSamples/CF)*SC
        y1_id8 = y2_id8 = DDDecBaselineEnd
        tl_green_wide.DrawLine(x1_id8,y1_id8,x2_id8,y2_id8)
        legend.AddEntry(tl_green_wide,"Baseline After Pulse","l")
        # Draw ellipses around the baseline lines since they are hard
        # to see sometimes when zoomed out.
        radius_x = abs(x1_id7-x1_id8)/10.0
        radius_y = abs(y1_id7-y1_id8)/10.0
        te_green.DrawEllipse(x1_id7, y1_id7, radius_x, radius_y, 0, 360,0)
        te_green.DrawEllipse(x2_id8, y1_id8, radius_x, radius_y, 0, 360,0)
    elif OutputPulseType in [8]: #Integral of deconvolved, baseline removed
        # These scaled Y-coordinates are used several times here.
        # They are just convenient places to start/end vertical lines.
        Y1_opt8, Y2_opt8 = -0.05*DDAmpl, 1.05*DDAmpl
        # id1source
        x1_id1 = x2_id1 = DD10pct*SC
        y1_id1, y2_id1 = Y1_opt8, Y2_opt8
        tl_green.DrawLine(x1_id1,y1_id1,x2_id1,y2_id1)
        legend.AddEntry(tl_green,"10% Amplitude","l")
        # id1bsource
        x1_id1b = x2_id1b = DD25pct*SC
        y1_id1b, y2_id1b = Y1_opt8, Y2_opt8
        tl_green_dashed.DrawLine(x1_id1b,y1_id1b,x2_id1b,y2_id1b)
        legend.AddEntry(tl_green_dashed,"25% Amplitude","l")
        # id2bsource
        x1_id2b = x2_id2b = DD75pct*SC
        y1_id2b, y2_id2b = Y1_opt8, Y2_opt8
        tl_green_dashed.DrawLine(x1_id2b,y1_id2b,x2_id2b,y2_id2b)
        legend.AddEntry(tl_green_dashed,"75% Amplitude","l")
        # id2source
        x1_id2 = x2_id2 = DD90pct*SC
        y1_id2, y2_id2 = Y1_opt8, Y2_opt8
        tl_green.DrawLine(x1_id2,y1_id2,x2_id2,y2_id2)
        legend.AddEntry(tl_green,"90% Amplitude","l")
        # id3source
        x1_id3 = x2_id3 = DDRealStopTime*SC
        y1_id3, y2_id3 = Y1_opt8, Y2_opt8
        tl_darkgray_dashed.DrawLine(x1_id3,y1_id3,x2_id3,y2_id3)
        legend.AddEntry(tl_darkgray_dashed,"Integration Stop Time","l")
        # id4source
        x1_id4 = x2_id4 = DDRealTriggerTime*SC
        y1_id4, y2_id4 = Y1_opt8, Y2_opt8
        tl_darkgray_dashed.DrawLine(x1_id4,y1_id4,x2_id4,y2_id4)
        legend.AddEntry(tl_darkgray_dashed,"Integration Start Time","l")
        # id7source
        x1_id7 = (DDRealTriggerTime - DecBaselineSamples/CF)*SC
        x2_id7 = DDRealTriggerTime*SC
        y1_id7 = y2_id7 = 0
        tl_green_wide.DrawLine(x1_id7,y1_id7,x2_id7,y2_id7)
        legend.AddEntry(tl_green_wide,"Baseline Before Pulse","l")
        # id8source
        x1_id8 = DDRealStopTime*SC
        x2_id8 = (DDRealStopTime + DecBaselineSamples/CF)*SC
        y1_id8 = y2_id8 = DDAmpl
        tl_green_wide.DrawLine(x1_id8,y1_id8,x2_id8,y2_id8)
        legend.AddEntry(tl_green_wide,"Baseline After Pulse","l")
        # Draw ellipses around the baseline lines since they are hard
        # to see sometimes when zoomed out.
        radius_x = abs(x1_id7-x1_id8)/10.0
        radius_y = abs(y1_id7-y1_id8)/10.0
        te_green.DrawEllipse(x1_id7, y1_id7, radius_x, radius_y, 0, 360,0)
        te_green.DrawEllipse(x2_id8, y1_id8, radius_x, radius_y, 0, 360,0)       
    elif OutputPulseType in [9]: 
        #Integral of deconvolved, baseline removed, in ADUs
        # These scaled Y-coordinates are used several times here.
        # They are just convenient places to start/end vertical lines.
        Y1_opt9, Y2_opt9 = -0.05*DDAmplADU, 1.05*DDAmplADU
        # id1source
        x1_id1 = x2_id1 = DD10pct*SC
        y1_id1, y2_id1 = Y1_opt9, Y2_opt9
        tl_green.DrawLine(x1_id1,y1_id1,x2_id1,y2_id1)
        legend.AddEntry(tl_green,"10% Amplitude","l")
        # id1bsource
        x1_id1b = x2_id1b = DD25pct*SC
        y1_id1b, y2_id1b = Y1_opt9, Y2_opt9
        tl_green_dashed.DrawLine(x1_id1b,y1_id1b,x2_id1b,y2_id1b)
        legend.AddEntry(tl_green_dashed,"25% Amplitude","l")
        # id2bsource
        x1_id2b = x2_id2b = DD75pct*SC
        y1_id2b, y2_id2b = Y1_opt9, Y2_opt9
        tl_green_dashed.DrawLine(x1_id2b,y1_id2b,x2_id2b,y2_id2b)
        legend.AddEntry(tl_green_dashed,"75% Amplitude","l")
        # id2source
        x1_id2 = x2_id2 = DD90pct*SC
        y1_id2, y2_id2 = Y1_opt9, Y2_opt9
        tl_green.DrawLine(x1_id2,y1_id2,x2_id2,y2_id2)
        legend.AddEntry(tl_green,"90% Amplitude","l")
        # id3source
        x1_id3 = x2_id3 = DDRealStopTime*SC
        y1_id3, y2_id3 = Y1_opt9, Y2_opt9
        tl_darkgray_dashed.DrawLine(x1_id3,y1_id3,x2_id3,y2_id3)
        legend.AddEntry(tl_darkgray_dashed,"Integration Stop Time","l")
        # id4source
        x1_id4 = x2_id4 = DDRealTriggerTime*SC
        y1_id4, y2_id4 = Y1_opt9, Y2_opt9
        tl_darkgray_dashed.DrawLine(x1_id4,y1_id4,x2_id4,y2_id4)
        legend.AddEntry(tl_darkgray_dashed,"Integration Start Time","l")
        # id7source
        x1_id7 = (DDRealTriggerTime - DecBaselineSamples/CF)*SC
        x2_id7 = DDRealTriggerTime*SC
        y1_id7 = y2_id7 = 0
        tl_green_wide.DrawLine(x1_id7,y1_id7,x2_id7,y2_id7)
        legend.AddEntry(tl_green_wide,"Baseline Before Pulse","l")
        # id8source
        x1_id8 = DDRealStopTime*SC
        x2_id8 = (DDRealStopTime + DecBaselineSamples/CF)*SC
        y1_id8 = y2_id8 = DDAmplADU
        tl_green_wide.DrawLine(x1_id8,y1_id8,x2_id8,y2_id8)
        legend.AddEntry(tl_green_wide,"Baseline After Pulse","l")
        # Draw ellipses around the baseline lines since they are hard
        # to see sometimes when zoomed out.
        radius_x = abs(x1_id7-x1_id8)/10.0
        radius_y = abs(y1_id7-y1_id8)/10.0
        te_green.DrawEllipse(x1_id7, y1_id7, radius_x, radius_y, 0, 360,0)
        te_green.DrawEllipse(x2_id8, y1_id8, radius_x, radius_y, 0, 360,0)
    # end the big if/elif block for OutputPulseType

    legend.Draw()
    
    ROOT.gPad.Modified()
    ROOT.gPad.Update()

    # Return the TGraph and the TLegend.  Other objects like the TLines
    # and TEllipses are leaked memory.
    return g,legend

# The following part is only used if you run the script as a program.

helpstring = """[Left/right arrow] keys go to the previous/next entry.
[Shift+arrow] skips 10  [Ctrl+Arrow] skips 100  [Ctrl+Shift+Arrow] skips 1000.
[Up/down arrow] keys change pulses within an event.
[u] switches between Samples and microseconds. 
[o] lets you select a new OutputPulseType. 
[j] lets you enter a new JSON input filename.
[z] fixes/unfixes the X-axis zoom level.
[g] jumps to a specific entry number.           [c] lets you enter a ROOT TCut.
[s] saves the canvas.   [r] redraws the figure. [q] quits.  
"""

outputhelp = """List of OutputPulseType:
[-1]: Raw pulse                           [0]: Baseline removed 
 [1]: Trapezoidal filter                  [2]: Baseline removed, light smoothing
 [3]: Deconvolved from Preamp             [4]: Deconvolved from Preamp, Smoothed
 [5]: Double Deconvolved, Harsh Smoothing [6], Double Deconvolved, Light Smoothing
 [7], Integral of Double Deconvolved Pulse
 [8], Baseline Removed Integral of Double Deconvolved Pulse
 [9], Baseline Removed Integral of Double Deconvolved Pulse (in ADUs)
 """
if __name__ == "__main__":
    import curses
    c1 = ROOT.TCanvas("c1","c1",800,600)

    def runloop(screen):
        global json_filename
        # Initialize some sane starting values before the loop.
        OutputPulseType = -1
        event_idx = 0
        entry_idx = 0 # Used for when a cut is entered.
        pulse_idx = 0
        xunits = ["Samples","microseconds"]
        xunits_i = 0
        fixed_xzoom = False
        cut = ""
        if T2: # If cut is "" the entrylist will have all the entries.
            T2_tree.Draw(">>entrylist",cut,"entrylist")
        else:
            T1_tree.Draw(">>entrylist",cut,"entrylist")
        entrylist = ROOT.gROOT.FindObject("entrylist")
        entrylistN = entrylist.GetN()
        screen.addstr(0,0,helpstring)
        screen.addstr(helpstring.count("\n"),0,outputhelp)
        statuslineidx = helpstring.count("\n") + outputhelp.count("\n")
        firstloop = True # Is this the first time through the event loop?
        while True:
            doupdate = True # Update the plot unless this changes.
            inc = 0 # Amount to increment the entry_idx.
            # Construct some status lines to print for the user.
            statusline1 = "Current OutputPulseType: %d, " % OutputPulseType +\
                         "event: %d, " % event_idx +\
                         "pulse: %d, " % pulse_idx + "entry: %d, " % entry_idx +\
                         "fixed_xzoom: %s" % str(fixed_xzoom)
            statusline2 = "JSON file: %s \n" % json_filename
            if T2:
                statusline3 = "T2 TCut: %s \n" % (cut if cut else "None")
            else:
                statusline3 = "T1 TCut: %s \n" % (cut if cut else "None")
            prompt_line = statuslineidx + 3
            error_line = prompt_line+1

            # Clear the status lines
            for i in range(4):
                screen.move(statuslineidx+i,0)
                screen.deleteln()
                screen.insertln()

            # Write the new status lines.
            screen.addstr(statuslineidx,0,statusline1)
            screen.addstr(statuslineidx+1,0,statusline2)
            screen.addstr(statuslineidx+2,0,statusline3)

            if firstloop:
                # If it's our first time through the loop, don't ask for 
                # user input, just draw the first entry.
                c = 'r'
                firstloop = False
            else:
                c = screen.getkey()

            # Now parse the user input.
            if c in ["q","Q"]: # Quit
                break
            elif c in ["u","U"]: # Change units.
                xunits_i = (xunits_i+1)%2 # Toggle between 0 and 1
            elif c in ["o","O"]: # Change OutputPulseType
                screen.addstr(prompt_line,0,"New OutputPulseType: ")
                curses.echo()
                newoutput = screen.getstr()
                curses.noecho()
                try: # Validate the number.
                    if not (-1 <= int(newoutput) <= 9):
                        raise ValueError
                    else:
                        OutputPulseType = int(newoutput)
                except:
                    errmsg = "OutputPulseType must enter an integer from -1 to 9.\n"
                    screen.addstr(error_line,0,errmsg)
            elif c in ['j',"J"]: # Change JSON file.
                screen.addstr(prompt_line,0,"Enter new JSON filename: ")
                curses.echo()
                newjson = screen.getstr()
                curses.noecho()
                if os.path.isfile(newjson):
                    json_filename = newjson
                else:
                    screen.addstr(error_line,0,"File does not exist.\n")
            elif c in ['g',"G"]: # Jump to a specific entry number.
                # Note that the number is going to be the entry in the
                # entrylist.  If a cut is applied, the entry number
                # does not correspond to the event number.
                screen.addstr(prompt_line,0,"Enter entry number: ")
                curses.echo()
                newevent = screen.getstr()
                curses.noecho()
                try:
                    entry_idx = int(newevent)
                except ValueError:
                    screen.addstr(error_line,0,"Entry numbers must be integers.\n")
            elif c in ['s','S']: # Save a figure.
                # Construct a file name.
                outname = "_".join(map(str,[runname,event_idx,pulse_idx,OutputPulseType]))
                fulloutname = os.path.join(output_dir,outname)
                with Quiet(): # This suppresses the ROOT save message.
                    ROOT.gPad.Modified()
                    ROOT.gPad.Update()
                    for of in output_formats:
                        c1.SaveAs(fulloutname+"."+of)
                msg = "Saved to %s.%s\n" % (fulloutname,str(output_formats))
                screen.addstr(error_line,0,msg)
                doupdate = False # Don't redraw just because it was saved.
            elif c in ['r','R']:
                # Do nothing, but this cycles the event loop.
                pass
            elif c in ['c','C']: # Enter a cut.
                if T2:
                    screen.addstr(prompt_line,0,"Enter T2 TCut: ")
                else:
                    screen.addstr(prompt_line,0,"Enter T1 TCut: ")
                curses.echo()
                cut = screen.getstr()
                curses.noecho()
                if T2:
                    T2_tree.Draw(">>entrylist",cut,"entrylist")
                else:
                    T1_tree.Draw(">>entrylist",cut,"entrylist")
                entrylist = ROOT.gROOT.FindObject("entrylist")
            elif c in ['z','Z']: # Fix/unfix zoomlevel
                if fixed_xzoom:
                    g.GetXaxis().UnZoom()
                    fixed_xzoom = False
                else:                    
                    xax = g.GetXaxis()
                    minx, maxx = xax.GetFirst(), xax.GetLast()
                    doupdate = False # Don't update from fixing zoom.
                    fixed_xzoom = True
            elif c == "KEY_LEFT": 
                inc = -1
            elif c == "KEY_RIGHT":
                inc = 1
            elif c == "KEY_SLEFT": # Shift + Arrow
                inc = -10
            elif c == "KEY_SRIGHT":
                inc = 10
            elif c == "kLFT5":     # Ctrl + Arrow
                inc = -100
            elif c == "kRIT5":
                inc = 100
            elif c == "kLFT6":     # Ctrl + Shift + Arrow
                inc = -1000
            elif c == "kRIT6":
                inc = 1000
            elif c == "KEY_UP":
                pulse_idx = min((pulse_idx + 1), num_pulses-1)
            elif c == "KEY_DOWN":
                pulse_idx = max((pulse_idx - 1), 0)
            else:
                screen.addstr(error_line,0,"Error, not a command: "+c)

            # Increment the event index but loop around nEntries.
            entry_idx = (entry_idx + inc)
            while entry_idx < 0:
                entry_idx += entrylistN
            entry_idx = entry_idx % entrylistN
            event_idx = entrylist.GetEntry(entry_idx)
            if doupdate:
                g,l = DDFilter(OutputPulseType,event_idx,pulse_idx,
                               json_filename,xunits[xunits_i])
                if fixed_xzoom:
                    g.GetXaxis().SetRange(minx,maxx)
                    ROOT.gPad.Modified()
                    ROOT.gPad.Update()


    # This runs the runloop function and handles curses setup/teardown.
    curses.wrapper(runloop)

##################################################################
# Description of output pulses and lines from analysis tutorial.
# OutputPulseType
# -1, "Raw pulse"
# 0, "Baseline Removed": 
    # the vertical green lines indicate the points where 50% of the
    # maximum is reached (the difference between the two is the width), 
    # the striped green line indicates the point where the maximum is 
    # reached, the vertical striped black lines indicate the points 
    # where the process thinks the pulse starts and ends.
# 1, "Trapezoidal Filter": 
    # the horizontal red line is the start threshold value, 
    # the vertical striped line is the resulting start point 
    # (note, we add a safety margin afterwards)
# 2, "Baseline Removed, Light Smoothing"
# 3, "Deconvolved from Preamp"
# 4, "Deconvolved from Preamp, Smoothed"
# 5, "Double Deconvolved, Harsh Smoothing": 
    # the horizontal black line is the end threshold, [de3 /de4 ->  
    # former is value set by user, latter is lowered value set by 
    # algo if former did not "trigger"]
    # the rightmost vertical striped line is the resulting end point [de2] 
    # (note, we add a safety margin afterwards)
# 6, "Double Deconvolved, Light Smoothing"
# 7, "Integral of Double Deconvolved Pulse"
# 8, "Baseline Removed Integral of Double Deconvolved Pulse": 
    # the vertical green solid (striped) lines indicate the point 
    # where the pulse reaches 10% and 90% [id1, id2] (25% and 75%) 
    # [id1b, id2b] of the total amplitude (the difference between 
    # the two is the rise time), 
    # the horizontal green lines indicates the value of the baseline 
    # before and after the pulse [id3, id4] (the difference between the two is
    # the amplitude).
# 9, "Baseline Removed Integral of Double Deconvolved Pulse (in ADUs)"

