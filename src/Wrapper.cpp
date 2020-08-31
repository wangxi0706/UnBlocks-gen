/*    
    UnBlocks-Gen: 3D rock mass generator and analyser
    Copyright (C) 2020 Leandro Lima Rasmussen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

void Generator_Wrapper();
void DFN_Wrapper();
void FractureSet_Wrapper();
void Fracture_Wrapper();
void Mapping_Wrapper();
void Block_Wrapper();

void UnBlocks_Greet(){
    std::cout <<  "UnBlocks-Gen: 3D rock mass generator and analyser" << std::endl;
    std::cout <<  "Copyright (C) 2020 Leandro Lima Rasmussen" << std::endl;
    std::cout <<  "This program comes with ABSOLUTELY NO WARRANTY." << std::endl;
    std::cout <<  "This is free software, and you are welcome to" << std::endl;
    std::cout <<  "redistribute it under certain conditions." << std::endl;
    std::cout <<  std::endl;                                                                                                                             
    std::cout << "                                                                                                                                 GGGGGGGGGGGGG                                       " << std::endl;  
    std::cout << "                                                                                                                               GGG::::::::::::G                                      " << std::endl;   
    std::cout << "                                                                                                                             GG:::::::::::::::G                                      " << std::endl;   
    std::cout << "                                                                                                                            G:::::GGGGGGGG::::G                                      " << std::endl;   
    std::cout << "                                                                                                                           G:::::G       GGGGGG    eeeeeeeeeeee    nnnn  nnnnnnnn    " << std::endl;    
    std::cout << "                                                                                                                          G:::::G                ee::::::::::::ee  n:::nn::::::::nn  " << std::endl;   
    std::cout << "                                                                                                                          G:::::G               e::::::eeeee:::::een::::::::::::::nn " << std::endl;   
    std::cout << "                                                                                                                          G:::::G    GGGGGGGGGGe::::::e     e:::::enn:::::::::::::::n" << std::endl;  
    std::cout << "                                                                                                                          G:::::G    G::::::::Ge:::::::eeeee::::::e  n:::::nnnn:::::n" << std::endl;  
    std::cout << "                                                                                                                          G:::::G    GGGGG::::Ge:::::::::::::::::e   n::::n    n::::n" << std::endl;  
    std::cout << "                                                                                                                          G:::::G        G::::Ge::::::eeeeeeeeeee    n::::n    n::::n" << std::endl;     
    std::cout << "                                                                                                                            G:::::GGGGGGGG::::Ge::::::::e            n::::n    n::::n" << std::endl;      
    std::cout << "                                                                                                                            G:::::GGGGGGGG::::Ge::::::::e            n::::n    n::::n" << std::endl;      
    std::cout << "UUUUUUUU     UUUUUUUU                 BBBBBBBBBBBBBBBBB   lllllll                                      kkkkkkkk              GG:::::::::::::::G e::::::::eeeeeeee    n::::n    n::::n" << std::endl;
    std::cout << "U::::::U     U::::::U                 B::::::::::::::::B  l:::::l                                      k::::::k                GGG::::::GGG:::G  ee:::::::::::::e    n::::n    n::::n" << std::endl;
    std::cout << "U::::::U     U::::::U                 B::::::BBBBBB:::::B l:::::l                                      k::::::k                   GGGGGG   GGGG    eeeeeeeeeeeeee    nnnnnn    nnnnnn" << std::endl;
    std::cout << "UU:::::U     U:::::UU                 BB:::::B     B:::::Bl:::::l                                      k::::::k                          " << std::endl;
    std::cout << " U:::::U     U:::::Unnnn  nnnnnnnn      B::::B     B:::::B l::::l    ooooooooooo       cccccccccccccccc k:::::k    kkkkkkk  ssssssssss   " << std::endl;
    std::cout << " U:::::D     D:::::Un:::nn::::::::nn    B::::B     B:::::B l::::l  oo:::::::::::oo   cc:::::::::::::::c k:::::k   k:::::k ss::::::::::s  " << std::endl;
    std::cout << " U:::::D     D:::::Un::::::::::::::nn   B::::BBBBBB:::::B  l::::l o:::::::::::::::o c:::::::::::::::::c k:::::k  k:::::kss:::::::::::::s " << std::endl;
    std::cout << " U:::::D     D:::::Unn:::::::::::::::n  B:::::::::::::BB   l::::l o:::::ooooo:::::oc:::::::cccccc:::::c k:::::k k:::::k s::::::ssss:::::s" << std::endl;
    std::cout << " U:::::D     D:::::U  n:::::nnnn:::::n  B::::BBBBBB:::::B  l::::l o::::o     o::::oc::::::c     ccccccc k::::::k:::::k   s:::::s  ssssss " << std::endl;
    std::cout << " U:::::D     D:::::U  n::::n    n::::n  B::::B     B:::::B l::::l o::::o     o::::oc:::::c              k:::::::::::k      s::::::s      " << std::endl;
    std::cout << " U:::::D     D:::::U  n::::n    n::::n  B::::B     B:::::B l::::l o::::o     o::::oc:::::c              k:::::::::::k         s::::::s   " << std::endl;
    std::cout << " U::::::U   U::::::U  n::::n    n::::n  B::::B     B:::::B l::::l o::::o     o::::oc::::::c     ccccccc k::::::k:::::k  ssssss   s:::::s " << std::endl;
    std::cout << " U:::::::UUU:::::::U  n::::n    n::::nBB:::::BBBBBB::::::Bl::::::lo:::::ooooo:::::oc:::::::cccccc:::::ck::::::k k:::::k s:::::ssss::::::s" << std::endl;
    std::cout << "  UU:::::::::::::UU   n::::n    n::::nB:::::::::::::::::B l::::::lo:::::::::::::::o c:::::::::::::::::ck::::::k  k:::::ks::::::::::::::s " << std::endl;
    std::cout << "    UU:::::::::UU     n::::n    n::::nB::::::::::::::::B  l::::::l oo:::::::::::oo   cc:::::::::::::::ck::::::k   k:::::ks:::::::::::ss  " << std::endl;
    std::cout << "      UUUUUUUUU       nnnnnn    nnnnnnBBBBBBBBBBBBBBBBB   llllllll   ooooooooooo       cccccccccccccccckkkkkkkk    kkkkkkksssssssssss    " << std::endl;
    std::cout <<  std::endl;                                                                                             
    std::cout << "Version 1.0" << std::endl; 
    std::cout << std::endl;
};

BOOST_PYTHON_MODULE(unblocks){
    UnBlocks_Greet();
    Generator_Wrapper();
    DFN_Wrapper();
    FractureSet_Wrapper();
    Fracture_Wrapper();
    Mapping_Wrapper();
    Block_Wrapper();
}
