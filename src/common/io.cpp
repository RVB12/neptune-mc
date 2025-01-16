
#include <fstream>

#include "eiron2/common/io.hpp"

namespace eiron2 {

    vector<index_t> load_intvector(std::string filename)
    {
        std::ifstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Could not open file " << filename << std::endl;
            exit(1);
        }

        index_t n;
        file >> n;

        vector<index_t> v(n);

        for (index_t i = 0; i < n; i++) {
            file >> v[i];
        }

        return v;
    }

    void save_intvector(vector<index_t> v, std::string filename)
    {
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Could not open file " << filename << std::endl;
            exit(1);
        }

        file << v.size() << std::endl;

        for (index_t i = 0; i < v.size(); i++) {
            file << v[i] << std::endl;
        }

        file.close();
    }


}