# pragma once

# include <iostream>
# include <vector>
# include <array>
# include <memory>
# include <stdexcept>

# include "OrdinateField.hpp"
# include "DiscreteOrdinateBC.hpp"
# include "PeriodicBC.hpp"
# include "DOMatrix.hpp"
# include "BoundaryCondition.hpp"
# include "PathConfig.hpp"

using namespace std ;

class OrdinatesConstructor{

    private:

        const SparseAddress& addr_ ;
        const FaceAddressedMesh3D& mesh_ ;
        const PathConfig& paths_ ;
        const WallProperty& wp_ ;
        const Field& temperatureField_ ;

        vector<unique_ptr<DOMatrix>>& intensity_ ;
        vector<OrdinateS> ordinateDirections_ ;
        int ordinate_ ;
        vector<OrdinateS> chosenOrdinate_ ;
        int n_ ;

        double sigma_ = 5.670374419e-8  ;
        double pi_ = 3.1415926535897932 ;

        // constant ordinates 
        vector<OrdinateS> s2 = {{0.5000000,0.7071068,0.5000000,1.5707963}} ;
        vector<OrdinateS> s4 = {{0.2958759,0.2958759,0.9082483,0.5235987},{0.2958759,0.9082483,0.2958759,0.5235987},{0.9082483,0.2958759,0.2958759,0.5235987}} ;
        vector<OrdinateS> s6 = {{0.1838670,0.1838670,0.9656013,0.1609517},{0.1838670,0.6950514,0.6950514,0.3626469},{0.1838670,0.9656013,0.1838670,0.1609517},{0.6950514,0.1838670,0.6950514,0.3626469},{0.6950514,0.6950514,0.1838670,0.3626469},{0.9656013,0.1838670,0.1838670,0.1609517}} ;
        vector<OrdinateS> s8 = {{0.1422555,0.1422555,0.9795543,0.1712359},{0.1422555,0.5773503,0.8040087,0.0992284},{0.1422555,0.8040087,0.5773503,0.0992284},{0.1422555,0.9795543,0.1422555,0.1712359},{0.5773503,0.1422555,0.8040087,0.0992284},{0.5773503,0.5773503,0.5773503,0.4617179},{0.5773503,0.8040087,0.1422555,0.0992284},{0.8040087,0.1422555,0.5773503,0.0992284},{0.8040087,0.5773503,0.1422555,0.0992284},{0.9795543,0.1422555,0.1422555,0.1712359}};

        // Field creator
        Field createIntensityField(DOMatrix& currentIntensity){

            SparseLinearVector internal_(mesh_.numCells()) ;
            vector<unique_ptr<BoundaryCondition>> boundary_ ;
            boundary_.resize(mesh_.boundaryPatches().size()) ;

            for(Patch p: mesh_.boundaryPatches()){

                const string temperatureBoundaryType = temperatureField_.boundary()[p.index()]->typeName();

                if(temperatureBoundaryType == "fixedvalue"){

                    map<int,double> faceVal_ ;
                    for(int faceIndex: p.faces()){
                        double val = wp_(p.index())*(sigma_*pow(temperatureField_.boundary()[p.index()]->boundaryValues()[faceIndex],4)/pi_) ;
                        faceVal_.insert({faceIndex, val}) ;
                    }
                    map<int,double>& wallTemperature = temperatureField_.boundary()[p.index()]->boundaryValues() ;
                    boundary_[p.index()] = make_unique<DiscreteOrdinateBC>(currentIntensity.I(), mesh_, currentIntensity.F(), wp_(p.index()), intensity_, faceVal_, wallTemperature) ;

                }
                else if(temperatureBoundaryType == "fixedgradient"){

                    map<int,double> faceVal_ ;
                    for(int faceIndex: p.faces()){
                        faceVal_.insert({faceIndex, 0.0}) ;
                    }
                    boundary_[p.index()] = make_unique<PeriodicBC>(currentIntensity.I(), mesh_, currentIntensity.F(), faceVal_) ;

                }
                else{
                    throw runtime_error("Unsupported temperature boundary type for DOM patch " + p.name() + ": " + temperatureBoundaryType);
                }
            }

            return Field(std::move(internal_), std::move(boundary_)) ;
        }

    public:

        OrdinatesConstructor(const SparseAddress& addr, const FaceAddressedMesh3D& mesh, PathConfig& paths, const WallProperty& wp, const Field& temperatureField, vector<unique_ptr<DOMatrix>>& intensity, int ordinate):addr_(addr),
                                        mesh_(mesh),
                                        paths_(paths),
                                        wp_(wp),
                                        temperatureField_(temperatureField),
                                        intensity_(intensity),
                                        ordinate_(ordinate)
        {
            createOrdinates() ;
            createOrdinateDirectionLinearSystems() ;
        }

        // non-const access to the DO matrices
        vector<unique_ptr<DOMatrix>>& intensity(){ return intensity_;}
        const vector<unique_ptr<DOMatrix>>& intensity() const{ return intensity_;}

        void createOrdinates(){
            switch(ordinate_){
                case 2:
                    n_ = s2.size() ;
                    chosenOrdinate_ = s2 ;
                    createOrdinateDirections() ;
                    break ;
                case 4:
                    n_ = s4.size() ;
                    chosenOrdinate_ = s4 ;
                    createOrdinateDirections() ;
                    break ;
                case 6:
                    n_ = s6.size() ;
                    chosenOrdinate_ = s6 ;
                    createOrdinateDirections() ;
                    break ;
                case 8:
                    n_ = s8.size() ;
                    chosenOrdinate_ = s8 ;
                    createOrdinateDirections() ;
                    break ;
                default:
                    throw std::invalid_argument("Unsupported DOM ordinate order S" + std::to_string(ordinate_) + ". Supported orders are S2, S4, S6 and S8.");
            }
        }

        void createOrdinateDirections(){

            for(int i=0; i<n_; ++i){

                OrdinateS arr = chosenOrdinate_[i] ;

                ordinateDirections_.push_back({arr.x,arr.y,arr.z,arr.w}) ;
                ordinateDirections_.push_back({-arr.x,arr.y,arr.z,arr.w}) ;
                ordinateDirections_.push_back({arr.x,-arr.y,arr.z,arr.w}) ;
                ordinateDirections_.push_back({arr.x,arr.y,-arr.z,arr.w}) ;
                ordinateDirections_.push_back({-arr.x,-arr.y,arr.z,arr.w}) ;
                ordinateDirections_.push_back({arr.x,-arr.y,-arr.z,arr.w}) ;
                ordinateDirections_.push_back({-arr.x,arr.y,-arr.z,arr.w}) ;
                ordinateDirections_.push_back({-arr.x,-arr.y,-arr.z,arr.w}) ;
                
            }
        }

        void createOrdinateDirectionLinearSystems(){

            for (std::size_t i = 0; i < ordinateDirections_.size(); ++i)
            {
                OrdinateS currentDirection = ordinateDirections_[i];
                auto currentMatrix = std::make_unique<DOMatrix>(addr_, mesh_, paths_, currentDirection, static_cast<int>(i));
                DOMatrix& currentIntensity = *currentMatrix;
                intensity_.push_back( std::move(currentMatrix) );
                Field intensityField = createIntensityField( currentIntensity );
                currentIntensity.getIntensityField( std::move(intensityField) );
                currentIntensity.I().clear();
                currentIntensity.I().initializeSolution(currentIntensity.iabf().internal());
            }
        }


        



} ;