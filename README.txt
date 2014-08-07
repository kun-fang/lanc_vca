About input files:



---------------
XXX.conf
---------------
This file is used to define all the sites of a cluster. It includes coordinates of cluster sites and
translational vectors. For example, the square.conf is following.

dimension of cluster:
2
number of sites:
4
coordinate
0.0	0.0
0.0	1.0
1.0	1.0
1.0	0.0
translate vector:
2.0	0.0
0.0	2.0

Here I include some most used clusters. If you need another cluster, just follow this format to create
your own cluster




----------------
cluster.input
----------------
This file include hopping information of a cluster and also a lattice. Sample.txt include some examples
of this file. for example, for a simple square-cluster-tiled square lattice:

square

lattice:
i j t Q1 Q2
1 2 1 0 0
1 2 1 0 -1
1 4 1 0 0
1 4 1 -1 0
2 3 1 0 0
2 3 1 -1 0
3 4 1 0 0
3 4 1 0 1

t: 1.00

mu: 2.00
U: 4.00
Tep: 0.00

cluster:
i j t
1 2 1
1 4 1
2 3 1
3 4 1

t: 1.00

mu: 2.00

SC:
delta: 1.0


The first line is the name of the cluster, make sure it is the same as the file name of the cluster
configuration file. After the line "lattice:" are lines that define the lattice and after the line 
"cluster:" are definitions about the cluster.

For hoppings of the lattice, there are five columns: i, j, t, Q1, Q2.
i, j are indices of cluster sites, make sure i<=j
t is the type of hopping, different types can be indicated below this part
Q1, Q2 are coefficents of translational vectors because sometimes the hopping will go out of the
cluster, the vector Q1*v1+Q2*v2 (v1,v2 are translational vectors defined in the XXX.conf file) will
lead the hopping to another cluster beside.

Following this part, you need to write all the hopping types line by line, like
t: 1.0
t': 2.0
t'': 3.0
the name of the hopping (t,t',t'') is up to you, there is a space between name and value

The next part is fixed for mu, U and temperature, don't change the order of them.

For hoppings of the cluster, there are only three columns: i, j, t
definations of them are the same as the lattice part. There is no translational vector for cluster.

The types of hoppings are for hoppings of cluster so can be different from the lattice

And then is mu value. This can also be different from the lattice. Do not remove this part

The last part is the part for Weiss fields that will be introduced into the cluster.
Write the name of the Weiss field first, for example: SC, AF
Then the symbol of the field and value of the field
The symbol is up to you. It is not important.
List the fields one by one if you need more than one field. However, 
                 ALWAYS PUT THE SC FIELD AT THE END

