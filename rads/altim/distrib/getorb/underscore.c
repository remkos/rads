/*
  Support routine for checking whether Fortran routines used in C programs
  or vise versa requires changing the C routine names, adding underscores
  or capitalizing the names.

  14-Feb-1998 - Remko Scharroo
*/

long underscore()
{
  return 0;
}

long UNDERSCORE()
{
  return 1;
}

long underscore_()
{
  return 2;
}

long _underscore()
{
  return 3;
}
