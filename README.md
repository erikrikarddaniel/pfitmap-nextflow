# pfitmap-nextflow
Nextflow pipeline for hierarchical protein classification


**profiles_hierarchy.tsv** file description:
A tsv file listing all hmm profiles that are used for hmm search described by the fields:
profile		Name of the hmm profile
prank		Rank of the profile (domain, family, class, subclass, group...)
psuperfamily	Name of the superfamily to which the profile belongs
pfamily		Name of the profile's family 
pclass		Name of the profile's class 
psubclass	Name of the profile's subclass 
pgroup		Name of the profile's group
psubgroup	Name of the profile's subgroup
version		Version of the profile's rank 
plen		Length of the profile (in bp)

The file must have a header describing the fields in the above order.

*example:*
```
profile	prank	psuperfamily	pfamily	pclass	psubclass	pgroup	psubgroup	version plen
NrdBh	subclass	Ferritin-like	NrdBR2lox	NrdB	NrdBh	0.6	337
```

