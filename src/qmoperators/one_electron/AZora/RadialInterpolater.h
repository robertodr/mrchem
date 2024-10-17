#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include "utils/PolyInterpolator.h"
#include <string>
#include <filesystem>

/**
 * @brief Read ZORA potential from file. Check if file exists and abort if it does not.
 * @param path Path to the file containing the ZORA potential
 * @param rGrid Vector containing the radial grid
 * @param vZora Vector containing the ZORA potential
 * @param kappa Vector containing the kappa parameter
*/
void readZoraPotential(const std::string path, Eigen::VectorXd &rGrid, Eigen::VectorXd &vZora, Eigen::VectorXd &kappa){
    std::vector<double> r, v, k;
    bool file_exists = std::filesystem::exists(path);
    if (!file_exists) {
        std::cerr << "File " << path << " does not exist." << std::endl;
        std::cout << "File " << path << " does not exist." << std::endl;
        exit(1);
    }
    std::ifstream file(path);
    std::string line;
    double r_, v_, k_;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        iss >> k_ >> r_ >> v_;
        r.push_back(r_);
        v.push_back(v_);
        k.push_back(k_);
    }
    file.close();
    rGrid = Eigen::Map<Eigen::VectorXd>(r.data(), r.size());
    vZora = Eigen::Map<Eigen::VectorXd>(v.data(), v.size());
    kappa = Eigen::Map<Eigen::VectorXd>(k.data(), k.size());
    // The kappa function is half of what is defined in the paper Scalar 
    // Relativistic Effects with Multiwavelets: Implementation and Benchmark
    // it is not used in the code, only the potential is used
}

class RadInterpolater {

    public:
    /**
     * @brief Construct a new Rad Interpolater object
     * @param element The element for which the ZORA potential is to be interpolated
     * @param data_dir The directory containing the ZORA potential data
     * @param mode The mode of interpolation. Either "potential" or "kappa"
    */
    RadInterpolater(const std::string element, std::string data_dir, const std::string mode){
        Eigen::VectorXd rGrid;
        Eigen::VectorXd vZora;
        Eigen::VectorXd kappa;

        this->mode = mode;
        std::string filename = data_dir + '/' + element + ".txt";

        readZoraPotential(filename, rGrid, vZora, kappa);
        if (mode == "kappa") {
            interpolator = std::make_shared<interpolation_utils::PolyInterpolator>(rGrid, kappa);
        } else if (mode == "potential") {
            interpolator = std::make_shared<interpolation_utils::PolyInterpolator>(rGrid, vZora);
        }

    }

    double evalf(const double &r) const {
        return interpolator->evalfLeftNoRightConstant(r);
    }

    protected:
    std::shared_ptr<interpolation_utils::PolyInterpolator> interpolator;
    std::string mode;

};