a = (e01 + e23);
b = (e02 + e31);
c = (e03 - e12);
c2 = (-e02 + e12 + e03);

apl = (e01 + e23) / sqrt(2);
bpl = (e02 + e31) / sqrt(2);
cpl = (e03 + e12) / sqrt(2);
am = (e01 - e23) / sqrt(2);
bm = (e02 - e31) / sqrt(2);
cm = (e03 - e12) / sqrt(2);

linea = (e02 + e12 + e03 - e31) / sqrt(2);
lineb = (-e01 + e12 + e03 + e23) / sqrt(2);
linec = (-e02 + e12 + e03 + e31) / sqrt(2);

S1 = e01 ^ e23;
S  = exp( (2 * S1) / 2);

T1 = e12 ^ e31;
T  = exp( T1 / 2 );
