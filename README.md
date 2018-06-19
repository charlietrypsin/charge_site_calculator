# charge_site_calculator
Function to create lists of protonation and deprotonation sites. For use with the Konermann charge hopping algorithm [1][2].

Takes a pdb file and outputs integer lists for use with pdb2gmx to assign charges to ionisable sites (LYS, ARG, ASP, GLU, HIS and both N- and C-termini) for a fully protonated state, fully deprotonated state and specific protonation state. Function arguments are filename and the specific protonation number.

1) Popa, V., Trecroce, D. A., McAllister, R. G., & Konermann, L. (2016). Collision-Induced Dissociation of Electrosprayed Protein Complexes: An All-Atom Molecular Dynamics Model with Mobile Protons. Journal of Physical Chemistry B, 120(23), 5114–5124. http://doi.org/10.1021/acs.jpcb.6b03035
2) Konermann, L. (2017). Molecular Dynamics Simulations on Gas-Phase Proteins with Mobile Protons: Inclusion of All-Atom Charge Solvation. Journal of Physical Chemistry B, 121(34), 8102–8112. http://doi.org/10.1021/acs.jpcb.7b05703
