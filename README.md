# ggff

A gff manipulation tool.

Basically I find that I regularly need to manipulate GFFs to clean up annotations or whatever.

There are a few parsers for GFFs but none that really handle the DAG structure particularly well.
Generally, the available GFF libraries will parse the file, but leave it up to the user how to deal with it.
My other project gffpal has done well for me for a while. Doing graph traversal is super easy and useful, but there are cases where it's more painful than i want.
Additionally, i wrote another parser as part of predector-utils that is actually slightly better, but isn't quite compatible with gffpal.

Here I'll start migrating to a proper graph backed structure using networkx.
I'll also support location based lookup using a lazily constructed interval tree.

The basic 9 column structure will be fixed, but the library will be generic over the attributes column.
Strict checking of the attributes according to the GFF3 spec can be great but it isn't useful if you need to clean up data or convert from GFF2 to GFF3 etc.

This will be a slow work in progress.
