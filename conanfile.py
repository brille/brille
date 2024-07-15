from conan import ConanFile
from conan.tools.cmake import cmake_layout


class BrilleRecipe(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeToolchain", "CMakeDeps"
    default_options = {
        "hdf5/*:shared": False,
        "hdf5/*:hl": False,
        "hdf5/*:with_zlib": False,
        "highfive/*:with_boost": False,
        "highfive/*:with_eigen": False,
        "highfive/*:with_xtensor": False,
        "highfive/*:with_opencv": False,
    }

    def requirements(self):
        self.requires("pybind11/2.13.1")
        # self.requires("hdf5/1.12.0")
        self.requires("catch2/3.6.0")
        self.requires("highfive/2.9.0")
        if self.settings.os == "Macos":
            self.requires("llvm-openmp/11.1.0")

    def layout(self):
        cmake_layout(self)
