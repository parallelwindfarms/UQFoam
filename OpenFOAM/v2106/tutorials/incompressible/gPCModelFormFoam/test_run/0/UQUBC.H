internalField   uniform (0 0 0);

boundaryField
{
    top
    {
        type            noSlip;
    }
    hills
    {
        type            noSlip;
    }
    "(inlet|outlet)"
    {
        type            cyclic;
    }
    "(front|back)"
    {
        type            empty;
    }
}
