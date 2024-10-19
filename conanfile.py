from conan import ConanFile
from conan.tools.files import copy
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps


class KFBaseRecipe(ConanFile):
    name = "kfbase"
    version = "1.0"
    generators = "CMakeToolchain"
    # Optional metadata
    license = "GPL-3.0"
    author = "Sergei Gribanov"
    url = "https://github.com/sergeigribanov/KFBase"
    description = "Base package for kinematic and vertex fitting"
    topics = ("kinematic fitting")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = "*"

    def requirements(self):
        self.requires("eigen/3.4.0")
        
    def deploy(self):
        copy(self, "*", src=self.package_folder, dst=self.deploy_folder)
    
    def layout(self):
        cmake_layout(self)

    # def generate(self):
    #     deps = CMakeDeps(self)
    #     deps.generate()
    #     tc = CMakeToolchain(self)
    #     tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["kfbase"]
