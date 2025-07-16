# Overview 
This is a table of [Clebsch-Gordan](https://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients) coefficients for SO(3)
and a [PARI/GP](https://pari.math.u-bordeaux.fr/) script to calculate them.

An equivalent and more comprehensive Java program is
in [arXiv:0908.3030](https://arxiv.org/src/0908.3030v4/anc/org/nevec/rjm/Wigner3j.java)

A Python program is in [arXiv:1102.5125](https://arxiv.org/abs/1102.5125)

# Invocation
To obtain PARI/GP use the Linux package manager, for example

```
zypper install pari-gp # openSUSE
sudo apt install pari-gp # Ubuntu
```

The ASCII table of the Clebsch-Gordan coefficients is then generated with

```
gp -q < CGord.gp
```

To increase or decrease the extend of the j<sub>1</sub> maximum quantum
number in the output, edit the last line of `CGord.gp`
and change the 9/2 to some other positive half-integer value.
Because PARI/GP does not have an interface to generic command line
arguments, this editing is apparently the only way to 
change that value.

The ASCII table of the 6j-Symbols is generated with

```
gp -q < 6j.gp
```

To increase or decrease the extend of the j<sub>1</sub> maximum quantum
number in the output, edit the last line of `6j.gp`
and change the 11/2 to some other positive half-integer value.

