.. _denovo_variant_generation:

Generating denovo variants
==========================
The actual logic of generating new variants is implemented in the plugins (Please see :ref:`how_to_write_variant_plugins`).
Once a list of proposed variants are received they are incorporated into the existing genome
using a zipper like operation (:py:func:`mitty.lib.variation.merge_variants`). This function is used often
throughout the code including for populations simulations (:ref:`pop_sims`).