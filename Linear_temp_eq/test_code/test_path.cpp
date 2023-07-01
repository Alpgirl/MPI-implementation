#include <mach-o/dyld.h>
#include <limits.h>

int main(int argc, char **argv)
{
  char buf [PATH_MAX];
  uint32_t bufsize = PATH_MAX;
  if(!_NSGetExecutablePath(buf, &bufsize))
    puts(buf);
  return 0;
}