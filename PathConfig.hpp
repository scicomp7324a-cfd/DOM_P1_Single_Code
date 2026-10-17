# pragma once 

#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <map>
#include <vector>
#include <string>

using namespace std ;
namespace fs = filesystem ;

class PathConfig{

    private:
        fs::path baseDir_ ;
        fs::path outputDir_ ;

        using MemFn = fs::path (PathConfig::*)() const;

    public:
        PathConfig()
            :baseDir_("Files"), outputDir_("../Output"){}

        PathConfig(const fs::path& baseDir, const fs::path& outputDir)
            :baseDir_(baseDir), outputDir_(outputDir){}

        // Base and Output directories 
        const fs::path& baseDir() const {return baseDir_;}
        const fs::path& outputDir() const {return outputDir_;}

        // Sub folders
        fs::path meshDir() const{return baseDir_ / "MeshFiles" ; } 
        fs::path fieldDir() const{return baseDir_ / "FieldFiles" ; }
        
        const fs::path getOutputDirectory(string parentFolder, string caseName, string subCaseName) const{
            fs::path mainDirectory = outputDir() / parentFolder / caseName / subCaseName ;
            return mainDirectory ;
        }

        const fs::path getOutputDirectory(string caseName, string subCaseName) const{
            fs::path mainDirectory = outputDir() / caseName / subCaseName ;
            return mainDirectory ;
        }

        const fs::path getOutputDirectory(string caseName) const{
            fs::path mainDirectory = outputDir() / caseName ;
            return mainDirectory ;
        }

        vector<string> getMeshFilenames(string caseName, string dimensions){
            fs::path meshPath = meshDir() / caseName / dimensions ;
            vector<string> filenames = {meshPath / "points.txt",
                                        meshPath / "faces.txt",
                                        meshPath / "cells.txt",
                                        meshPath / "boundary.txt",
                                        meshPath / "dimensions.txt"} ;
            return filenames ;
        }

        vector<string> getMeshFilenames(string parentFolder, string caseName, string dimensions){
            fs::path meshPath = meshDir() / parentFolder / caseName / dimensions ;
            vector<string> filenames = {meshPath / "points.txt",
                                        meshPath / "faces.txt",
                                        meshPath / "cells.txt",
                                        meshPath / "boundary.txt",
                                        meshPath / "dimensions.txt"} ;
            return filenames ;
        }

        // Field files based on name 
        fs::path UxFieldFile() const{ return fieldDir() / "PressureVelocity" / "Ux.txt" ;}
        fs::path UyFieldFile() const{ return fieldDir() / "PressureVelocity" / "Uy.txt" ;}
        fs::path UzFieldFile() const{ return fieldDir() / "PressureVelocity" / "Uz.txt" ;}
        fs::path PFieldFile() const{ return fieldDir() / "PressureVelocity" / "P.txt" ;}

        fs::path TFieldFile() const{ return fieldDir() / "T.txt" ;}
        fs::path GFieldFile() const{ return fieldDir() / "G.txt" ;}

        fs::path getTFieldFile(string parentFolder, string caseName, string subCaseName){
            return fieldDir() / parentFolder / caseName / subCaseName / "T.txt" ;
        }

        fs::path getGFieldFile(string parentFolder, string caseName, string subCaseName){
            return fieldDir() / parentFolder / caseName / subCaseName / "G.txt" ;
        }

        fs::path getTFieldFile(string caseName, string subCaseName){
            return fieldDir() / caseName / subCaseName / "T.txt" ;
        }

        fs::path getGFieldFile(string caseName, string subCaseName){
            return fieldDir() / caseName / subCaseName / "G.txt" ;
        }

        fs::path getTFieldFile(string caseName){
            return fieldDir() / caseName / "T.txt" ;
        }

        fs::path getGFieldFile(string caseName){
            return fieldDir() / caseName / "G.txt" ;
        }

        void validate() const {
            if (!fs::exists(baseDir_)) {
                throw std::runtime_error("Base directory not found: " + baseDir_.string());
            }
            if (!fs::exists(meshDir())) {
                throw std::runtime_error("Mesh directory not found: " + meshDir().string());
            }
            if (!fs::exists(fieldDir())) {
                throw std::runtime_error("Field directory not found: " + fieldDir().string());
            }
        }

        void ensureOutputDir() const {
            fs::create_directories(outputDir_);
        }

        void createOutputDirectory(string caseName) const {
            fs::path dir = (*this).getOutputDirectory(caseName) ;
            fs::create_directories(dir) ;
        }

        void createOutputDirectory(string caseName, string subCaseName) const {
            fs::path dir = (*this).getOutputDirectory(caseName, subCaseName) ;
            fs::create_directories(dir) ;
        }

        void createOutputDirectory(string parentFolder, string caseName, string subCaseName) const {
            fs::path dir = (*this).getOutputDirectory(parentFolder, caseName, subCaseName) ;
            fs::create_directories(dir) ;
        }

} ;

