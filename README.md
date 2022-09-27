# PlasmidActionIBM
This is a model about plasmid dynamics in surface-attached isogenic bacterial populations. It can be used to explore the faith of a single or of two interacting, very distinct types of plasmids: non-transmissible and transmissible (conjugative) plasmids; and shows how transmissible plasmids indirectly help non-transmissible plasmids to persist.

Hosts with non-transmissible plasmids are able to outcompete hosts with transmissible plasmids (as the transmissible plasmids are more costly and more easily lost from a host with both plasmids). Hosts with transmissible plasmids can lead to decreases in the density of plasmid-free cells (through conjugation). Finally, plasmid-free cells outcompete hosts with non-transmissible plasmids (due to the cost of these plasmids). The upshot of this non-transitive dynamic is the maintenance of the non-transmissible plasmid.

In addition to unraveling this very general mechanism for plasmid coexistence, the model demonstrates how local conditions triggering switches in the individual transfer competence of bacteria bearing transmissible plasmids affect population dynamics. By accounting for such quorum sensing mechanisms, the model can be considered more realistic than most models simulating conjugation dynamics and is believed to be able to serve as a very sophisticated model platform for many additonal future studies on this topic.

The model is implemented in the multi-agent programmable modeling environment NetLogo, which is freely available from https://ccl.northwestern.edu/netlogo/. It was created under NetLogo Version 6.0 (submitted as appendix to the article below) and refurbished in NetLogo 6.2.2, including a ‘reset’-procedure to set the model parameters to default values and including a video option using the vid-extension.

Downloading and running the model file with NetLogo allows you to reproduce some predefinded simulation experiments, stored as setups in the 'BehaviorSpace' as well as to run simulations with self-defined parameters.

A complete description of the model (ODD protocol) can be found in this article introducing the model:
Werisch, M., Berger, U. Berendonk, T. (2017): Conjugative plasmids enable the maintenance of low cost non-transmissible plasmids. Plasmid 91: 96-104. https://doi.org/doi:10.1016/j.plasmid.2017.04.004

Please indicate any use of this model that contributes to a publication with a reference to this article. If you plan to use or modify this model for your own research, please contact martin.zwanzig@tu-dresden.de to avoid working on the same or very similar research questions in parallel.
