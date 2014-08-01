SIMD_HartreeFock

This is a very rough not working start of an interface to the 
Libint integral engine. From here it will eventaully turn into a 
fully vectorized direct Hartree Fock code.

The untilmate goal is a direct Fock built 
that takes advantage of the Vectorized integral engine.
As it sits it does not work and debugging features have made it 
painfully slow. Once it is working those will be removed and 
the perormance will be restored.
Libint and eigen are required to run this. There exists no make file

Jonathan Dullea
jonathan.dullea@gmail.com
