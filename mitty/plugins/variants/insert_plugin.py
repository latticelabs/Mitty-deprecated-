"""This is the stock insertion generator.

Note: This never generates a deletion at the first base of a sequence.

"""
__example_param_text = """
{
  "chromosome": [1],  # Which chromosomes to run the variant generator on
  "phet": 0.8,        # Probability of heterozygous variation
  "p": 0.0001,        # Per-base probability of having an insertion
  "t_mat": [[ 0.32654629,  0.17292732,  0.24524503,  0.25528135],  # Base transition matrix
            [ 0.3489394,   0.25942695,  0.04942584,  0.3422078],
            [ 0.28778188,  0.21087004,  0.25963262,  0.24171546],
            [ 0.21644706,  0.20588717,  0.24978216,  0.32788362]],
  "p_end": 0.1,       # Probability of chain ending
  "max_len": 1000     # Maximum length of insertion
}
"""

_description = """
Stock insertion model that generates sequences with same base transition matrix as the human genome and creates a
power-law distribution of insertion lengths.
A typical parameter set resembles
""" + __example_param_text

_example_params = eval(__example_param_text)


from mitty.lib.variation import new_variation
import mitty.lib.util as mutil
import logging
logger = logging.getLogger(__name__)


def variant_generator(ref={},
                      chromosome=None,
                      p=0.01,
                      phet=0.5,
                      t_mat=None,
                      p_end=0.1,
                      max_len=1000,
                      master_seed=1,
                      **kwargs):
  assert 0 <= p <= 1.0, "Probability out of range"
  assert 0 <= phet <= 1.0, "Probability out of range"
  assert 0 <= p_end <= 1.0, "Probability out of range"
  if t_mat is None:
    t_mat = [[0.32654629, 0.17292732, 0.24524503, 0.25528135],
             [0.3489394, 0.25942695, 0.04942584, 0.3422078],
             [0.28778188, 0.21087004, 0.25963262, 0.24171546],
             [0.21644706, 0.20588717, 0.24978216, 0.32788362]]

  logger.debug('Master seed: {:d}'.format(master_seed))
  ins_loc_rng, ins_markov_rng, het_rng, copy_rng = mutil.initialize_rngs(master_seed, 4)

  pt_mat = mutil.add_p_end_to_t_mat(t_mat, p_end)

  for chrom in chromosome:
    ref_seq = ref[chrom]
    #ins_locs, = numpy.nonzero(ins_loc_rng.rand(len(ref_seq)) < p)
    ins_locs = [x for x in mutil.place_poisson(ins_loc_rng, p, len(ref_seq)) if ref_seq[x] != 'N']
    if len(ins_locs) == 0:
      continue  # No variants here
    ins_list = mutil.markov_sequences(ref_seq, ins_locs, max_len, pt_mat, ins_markov_rng)
    het_type = mutil.zygosity(len(ins_locs), phet, het_rng, copy_rng)

    yield {chrom: [new_variation(pos + 1, pos + 2, ref_seq[pos], ref_seq[pos] + the_insert, het)
           for het, pos, the_insert in zip(het_type, ins_locs, ins_list)]}


def test():
  """Should produce variants"""
  from mitty.lib.variation import Variation
  ref = {
    1: 'ACTGACTGACTG',
    2: 'ACTGACTGACTGACTGACTGACTG'
  }
  vg = variant_generator(ref, chromosome=[1, 2],
                         t_mat=[[0.32654629, 0.17292732, 0.24524503, 0.25528135],  # Base transition matrix
                                [0.3489394, 0.25942695, 0.04942584, 0.3422078],
                                [0.28778188, 0.21087004, 0.25963262, 0.24171546],
                                [0.21644706, 0.20588717, 0.24978216, 0.32788362]])
  for v_list in vg:
    assert isinstance(v_list.values()[0][0], Variation), v_list


if __name__ == "__main__":
    print _description