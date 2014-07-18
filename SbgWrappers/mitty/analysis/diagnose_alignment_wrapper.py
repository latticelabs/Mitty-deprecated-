"""This is the wrapper for the analysis program diagnose_alignment.py

Command line parameters are

plot_bad_alignments.py --file=F
"""
from sbgsdk import define, Process, require
import os


#@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class PlotAlign(define.Wrapper):
  class Inputs(define.Inputs):
    data_file = define.input(name='Data file (.pkl)', description='.pkl file of table of read analysis')

  class Outputs(define.Outputs):
    outfig = define.output(name='Figure file', description='Summary figure of alignment diagnostics')

  def execute(self):
    p = Process('python', '/Mitty/Analysis/plot_bad_alignments.py', '--file', self.inputs.data_file)
    p.run()

    self.outputs.outfig = os.path.splitext(self.inputs.data_file)[0] + '_analyzed.pdf'
    self.outputs.outfig.meta = self.inputs.data_file.make_metadata(file_type='.pdf')