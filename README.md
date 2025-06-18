# Overview 
This is a table of Clebsch-Gordan coefficients for SO(3)
and a [PARI/GP](https://pari.math.u-bordeaux.fr/) script to calculate them.

# Invocation
To obtain PARI/GP use the local Linux package manager, for example

```
zypper install pari-gp # openSUSE
sudo apt install pari-gp # Ubuntu
```

The ASCII table of the Clebsch-Gordan coefficients is then generated with

```
gp -q < CGord.gp
```

To increase or decrease the extend of the j1 maximum quantum
number in the output, edit the last line of `CGord.gp`
and change the 9/2 to some other positive half-integer value.
Because PARI/GP does not have an interface to generic command line
arguments, this editing is apparently the only mean to 
change that value.

The ASCII table of the 6j-Symbols is then generated with

```
gp -q < 6j.gp
```

To increase or decrease the extend of the j1 maximum quantum
number in the output, edit the last line of `6j.gp`
and change the 11/2 to some other positive half-integer value.

