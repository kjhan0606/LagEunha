#include <quadmath.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc,char**argv){
    try{
        if(argc<2)throw std::runtime_error("at least one sealed wire CSV required");
        size_t rows=0,values=0;
        for(int file=1;file<argc;file++){
            std::ifstream in(argv[file]);
            if(!in)throw std::runtime_error("read wire CSV");
            std::string line;
            if(!std::getline(in,line)||line.find("event_lower,event_upper")==std::string::npos)
                throw std::runtime_error("expected source-edge wire header");
            while(std::getline(in,line)){
                std::vector<std::string>field;
                size_t first=0;
                for(size_t next;(next=line.find(',',first))!=std::string::npos;first=next+1)
                    field.push_back(line.substr(first,next-first));
                field.push_back(line.substr(first));
                if(field.size()!=22)throw std::runtime_error("complete wire CSV row");
                for(int index: {11,12}){
                    const auto&word=field[index];
                    char*end=nullptr;
                    __float128 value=strtoflt128(word.c_str(),&end);
                    if(end!=word.c_str()+word.size()||!finiteq(value))
                        throw std::runtime_error("finite binary128 event endpoint");
                    char output[128];
                    int printed=quadmath_snprintf(output,sizeof(output),"%.36Qg",value);
                    if(printed<=0||size_t(printed)>=sizeof(output)||word!=output)
                        throw std::runtime_error("binary128 event decimal did not roundtrip");
                    values++;
                }
                rows++;
            }
            if(!in.eof())throw std::runtime_error("wire CSV read error");
        }
        std::cout<<"BINARY128_WIRE_ROUNDTRIP_PASS rows="<<rows
            <<" endpoints="<<values<<std::endl;
        return 0;
    }catch(const std::exception&e){
        std::cerr<<"BINARY128_WIRE_ROUNDTRIP_FAIL "<<e.what()<<std::endl;
        return 2;
    }
}
