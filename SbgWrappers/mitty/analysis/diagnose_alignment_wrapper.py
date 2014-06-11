"""This is the wrapper for the analysis program diagnose_alignment.py

Command line parameters are

diagnose_alignment.py --correctbam=CBAM --wrongbam=WBAM  --smalla=SMALLA  --outpkl=OUTPKL  --outfig=OUTFIG
"""
from sbgsdk import define, Process, require
import os


@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class PlotAlign(define.Wrapper):
  class Inputs(define.Inputs):
    good_bam = define.input(name='good BAM', description='.bam (and index) file of correctly aligned reads')
    bad_bam = define.input(name='bad BAM', description='.bam (and index) file of incorrectly aligned reads')
    smalla = define.input(name='SMALLA', description='The reference sequence .smalla file', required=True)

  class Outputs(define.Outputs):
    outpkl = define.output(name='STATS file', description='.pkl file of alignment diagnostics')
    outfig = define.output(name='Figure file', description='Summary figure of alignment diagnostics')

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    correct_bam = self.inputs.good_bam
    wrong_bam = self.inputs.bad_bam
    smallafname = self.inputs.smalla

    outpkl = os.path.join(output_dir, os.path.basename(correct_bam).replace('_correct.bam', '_alignment_diagnostic.pkl'))
    outfig = os.path.join(output_dir, os.path.basename(correct_bam).replace('_correct.bam', '_alignment_diagnostic_fig.png'))

    p = Process('python', '/Mitty/Analysis/diagnose_alignment.py',
                '--correctbam', correct_bam,
                '--wrongbam', wrong_bam,
                '--smalla', smallafname,
                '--outpkl', outpkl,
                '--outfig', outfig)
    p.run()

    self.outputs.outpkl = outpkl
    self.outputs.outpkl.meta = self.inputs.smalla.make_metadata(file_type='Python pickle')

    self.outputs.outfig = outfig
    self.outputs.outfig.meta = self.inputs.smalla.make_metadata(file_type='png')