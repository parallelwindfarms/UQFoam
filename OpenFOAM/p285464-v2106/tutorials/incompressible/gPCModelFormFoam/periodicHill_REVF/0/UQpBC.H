internalField   uniform 0;

boundaryField
{
    top
    {
        type            zeroGradient;
    }
    hills
    {
        type            zeroGradient;
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


