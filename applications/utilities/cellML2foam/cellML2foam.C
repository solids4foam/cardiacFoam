/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cardiacFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
  Thin wrapper to launch the cellML2foam Python CLI
\*---------------------------------------------------------------------------*/

#include <filesystem>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

int main(int argc, char** argv)
{
    namespace fs = std::filesystem;

    // Resolve script path relative to this source file location
    fs::path srcFile(__FILE__);
    fs::path script = srcFile.parent_path() / "scripts" / "cli.py";

    if (!fs::exists(script))
    {
        std::cerr << "cellML2foam: script not found: " << script << "\n";
        return 1;
    }

    std::vector<char*> args;
    args.reserve(static_cast<size_t>(argc) + 2);
    args.push_back(const_cast<char*>("python"));
    args.push_back(const_cast<char*>(script.c_str()));
    for (int i = 1; i < argc; ++i)
    {
        args.push_back(argv[i]);
    }
    args.push_back(nullptr);

    execvp("python", args.data());
    std::cerr << "cellML2foam: failed to exec python\n";
    return 1;
}
