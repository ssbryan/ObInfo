// instantiate an ObjInput object with cmdline filename
// that's it

#include "ObjInput.h"

#include <string>

int main(int argc, char* argv[])
{
    // read in initialization file, create Objects
    std::string infname;
    std::string outfname;

    // get init file from args (default "init.dat")
    if (argc == 2)
    {
        infname = argv[1];
    }
    
    if ((argc != 2) || !infname.compare("-?") || !infname.compare("/?") || !infname.compare("-help") || !infname.compare("/help"))
    {
        printf("Usage:  obinfo inputfilename\n");
        return 1;
    }
	
    outfname = infname;
    size_t n = outfname.rfind('.');

    if (n != std::string::npos)
    {
        size_t diff = outfname.length() - n;
        outfname.replace (n, diff, ".info");

        if (outfname.compare(infname) == 0)
        {
            // yell at user and tell them to use  different input extension
            printf("Input extension:  Do not use \".ascii\" since that's the default output extension\n");
            return 1;
        }
    }
    else
    {
        outfname += std::string(".ascii");
    }

    ObjInput oinp(infname, outfname);
    return oinp.ReadAndWriteInfo();
}
