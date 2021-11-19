//
// Created by g on 19/11/2021.
//

#ifndef BRILLE_HDF_INTERFACE_HPP
#define BRILLE_HDF_INTERFACE_HPP

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace brille {

    template<class T>
    bool lists_to_hdf(const std::vector<std::vector<T>>& lists, const std::string& filename, const std::string& entry){
        using namespace HighFive;
        File file(filename, File::ReadWrite | File::Create | File::Truncate);
        for (size_t i=0; i<lists.size(); ++i){
            file.createDataSet<T>(entry+"/"+std::to_string(i), DataSpace::From(lists[i]));
        }
        // write-out the length just to be safe?
        return true;
    }

    template<class T>
    std::vector<std::vector<T>> lists_from_hdf(const std::string& filename, const std::string& entry, size_t no){
        using namespace HighFive;
        File file(filename, File::ReadOnly);
        std::vector<std::vector<T>> out(no);
        for (size_t i=0; i<no; ++i){
            DataSet dataset = file.getDataSet(entry + "/" + std::to_string(i));
            dataset.read(out[i]);
        }
        return out;
    }

}
#endif //BRILLE_HDF_INTERFACE_HPP
