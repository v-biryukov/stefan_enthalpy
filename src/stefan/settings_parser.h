#pragma once
#include <stdlib.h> 
#include <iostream>
#include <sstream>
#include <array>
#include <vector>
#include <string>
#include <algorithm>
#include "../../3rdparty/tinyxml/tinyxml.h"
#include "../../3rdparty/tinyxml/tinystr.h"

struct Settings
{
    int dims_count;
    struct MeshSettings
    {
        std::string mesh_file;

        struct MediumParams
        {
            double T_phase, density_L, density_S;
            double thermal_conductivity_L, thermal_conductivity_S;
            double specific_heat_fusion, specific_heat_capacity_L, specific_heat_capacity_S;
        };
        std::vector<MediumParams> medium_params_info;
        
        struct BoundarySettings
        {
            enum BoundarySettingTypes
            {
                Custom,
                FixedFlux,
                FixedTemperature
            };
            BoundarySettingTypes type;
            int axis;
            int side;
            std::array<double, 3> params;
        };
        std::vector<BoundarySettings> boundary_settings_info;
    } mesh_settings;

    struct SnapshotSettings
    {
        int period_frames;
        std::string snaphot_data_path;
        std::vector<int> stride;
        bool write_temperature;
        bool write_enthalpy;
        bool write_thermal_elasticity;
        bool write_state;
    };
    SnapshotSettings snapshot_settings_info;

    struct TaskSettings
    {
        int number_of_steps;
        double time_step;
        struct InitialStateSettings
        {
            enum InitialStateSettingsType
            {
                PerNode,
                PerSubmesh
            };
            InitialStateSettingsType type;
            std::string initial_state_data_file;
            std::vector<double> initial_temperatures_by_submesh;
        } initial_state_settings;
    } task_settings;
};


Settings ParseFile(std::string xml_file_path);
void ParseMeshInfo(TiXmlElement* mesh_info_element, Settings& settings);
void ParseSnapshotInfo(TiXmlElement* snapshot_info_element, Settings& settings);
void ParseTaskInfo(TiXmlElement* task_info_element, Settings& settings);

template <typename T>
int ParseVector(TiXmlElement* element, const std::string name, std::vector<T>* vector)
{
    std::string text_value;
    text_value = std::string(element->Attribute(name.c_str()));
    std::stringstream iss( text_value );
    (*vector).clear();
    T temp_number;
    while ( iss >> temp_number )
        (*vector).push_back( temp_number );
    return strspn( text_value.c_str(), " -.0123456789" ) == text_value.size();
}

int ParseScalar(TiXmlElement* element, const std::string name, float* value)
{
    return element->QueryFloatAttribute(name.c_str(), value);
}

int ParseScalar(TiXmlElement* element, const std::string name, double* value)
{
    return element->QueryDoubleAttribute(name.c_str(), value);
}

int ParseBool(TiXmlElement* element, const std::string name, bool* value)
{
    return element->QueryBoolAttribute(name.c_str(), value);
}

void ParseString(TiXmlElement* element, const std::string name, std::string* value)
{
    *value = std::string(element->Attribute(name.c_str()));
}

template<typename T>
int ParseUnsigned(TiXmlElement* element, const std::string name, T* value)
{
    int signedValue = *value;
    int res = element->QueryIntAttribute(name.c_str(), &signedValue);
    *value = signedValue;
    return res;
}


Settings ParseFile(std::string xml_file_path)
{
    Settings settings;
    TiXmlDocument xml_doc;

    if (!xml_doc.LoadFile(xml_file_path.c_str()))
    {
        std::cerr << "Error: Loading " << xml_file_path << " fails with following error: " <<
            std::string(xml_doc.ErrorDesc()) << " in row " << xml_doc.ErrorRow() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    TiXmlElement* xml_element = xml_doc.FirstChildElement("Settings");
    if (!xml_element)
    {
        std::cerr << "Error: There is no Settings element \n";
        std::exit(EXIT_FAILURE);
    }

    if(xml_element->QueryIntAttribute("dimsCount", &settings.dims_count) != TIXML_SUCCESS)
    {
        std::cerr << "Error: There is no dimsCount element \n";
        std::exit(EXIT_FAILURE);
    }

    TiXmlElement* mesh_info_element = xml_element->FirstChildElement("Mesh");
    if (mesh_info_element)
    {
        ParseMeshInfo(mesh_info_element, settings);
    } else
    {
        std::cerr << "Error: Error: There is no Mesh section \n";
        std::exit(EXIT_FAILURE);
    }

    TiXmlElement* snapshot_info_element = xml_element->FirstChildElement("Snapshot");
    if (snapshot_info_element)
    {
        ParseSnapshotInfo(snapshot_info_element, settings);
    } else
    {
        std::cerr << "Error: There is no Snapshot section \n";
        std::exit(EXIT_FAILURE);
    }


    TiXmlElement* task_info_element = xml_element->FirstChildElement("Task");
    if (task_info_element)
    {
        ParseTaskInfo(task_info_element, settings);
    } else
    {
        std::cerr << "Error: There is no Task section \n";
        std::exit(EXIT_FAILURE);
    }
    return settings;

}

void ParseBoundaryAxis(TiXmlElement* local_boundaries_element, int& axis)
{
    std::string axis_temp_string;
    ParseString(local_boundaries_element, "axis", &axis_temp_string);
    if (!axis_temp_string.compare("x"))
            axis = 0;
    else if (!axis_temp_string.compare("y"))
            axis = 1;
    else if (!axis_temp_string.compare("z"))
            axis = 2;
    else
    {
            std::cerr << "Error: Wrong boundary axis settings \n";
            std::exit(EXIT_FAILURE);
    }
}

void ParseMeshInfo(TiXmlElement* mesh_info_element, Settings& settings)
{
    ParseString(mesh_info_element, "meshFile", &settings.mesh_settings.mesh_file);
    TiXmlElement* medium_params_element = mesh_info_element->FirstChildElement("MediumParams");
    if (medium_params_element)
    {

        TiXmlElement * submesh_element = medium_params_element->FirstChildElement("Submesh");
        while (submesh_element)
        {
            int index;
            Settings::MeshSettings::MediumParams temp_medium_params;
            ParseScalar(submesh_element, "TPhase", &(temp_medium_params.T_phase));
            ParseScalar(submesh_element, "densityL", &(temp_medium_params.density_L));
            ParseScalar(submesh_element, "densityS", &(temp_medium_params.density_S));
            ParseScalar(submesh_element, "thermalConductivityL", &(temp_medium_params.thermal_conductivity_L));
            ParseScalar(submesh_element, "thermalConductivityS", &(temp_medium_params.thermal_conductivity_S));
            ParseScalar(submesh_element, "specificHeatFusion", &(temp_medium_params.specific_heat_fusion));
            ParseScalar(submesh_element, "specificHeatCapacityL", &(temp_medium_params.specific_heat_capacity_L));
            ParseScalar(submesh_element, "specificHeatCapacityS", &(temp_medium_params.specific_heat_capacity_S));
            settings.mesh_settings.medium_params_info.push_back(temp_medium_params);
            submesh_element = submesh_element->NextSiblingElement("Submesh");
        };
    }

    settings.mesh_settings.boundary_settings_info.resize(2 * settings.dims_count);
    // Setting default boundary conditions: fixed flux = 0.0
    for (int i = 0; i < settings.mesh_settings.boundary_settings_info.size(); ++i)
    {
        settings.mesh_settings.boundary_settings_info[i].params = {0, 1, 0};
    }
    // Reading boundary conditions settings
    TiXmlElement* boundaries_element = mesh_info_element->FirstChildElement("Boundaries");
    if (boundaries_element)
    {
        TiXmlElement * custom_boundaries_element = boundaries_element->FirstChildElement("Custom");
        while (custom_boundaries_element)
        {
                Settings::MeshSettings::BoundarySettings temp_boundary_setings;
                temp_boundary_setings.type = Settings::MeshSettings::BoundarySettings::BoundarySettingTypes::Custom;
                ParseBoundaryAxis(custom_boundaries_element, temp_boundary_setings.axis);
                ParseUnsigned(custom_boundaries_element, "side", &temp_boundary_setings.side);
                std::vector<double> temp_vector;
                ParseVector(custom_boundaries_element, "params", &temp_vector);
                std::copy_n(temp_vector.begin(), 3, temp_boundary_setings.params.begin());
                settings.mesh_settings.boundary_settings_info[2*temp_boundary_setings.axis + temp_boundary_setings.side] = temp_boundary_setings;
                custom_boundaries_element = custom_boundaries_element->NextSiblingElement("Custom");
        }

        TiXmlElement* fixedflux_boundaries_element = boundaries_element->FirstChildElement("FixedFlux");
        while (fixedflux_boundaries_element)
        {
            Settings::MeshSettings::BoundarySettings temp_boundary_setings;
            temp_boundary_setings.type = Settings::MeshSettings::BoundarySettings::BoundarySettingTypes::FixedFlux;
            ParseBoundaryAxis(fixedflux_boundaries_element, temp_boundary_setings.axis);
            ParseUnsigned(fixedflux_boundaries_element, "side", &temp_boundary_setings.side);
            double flux;
            ParseScalar(fixedflux_boundaries_element, "flux", &(flux));
            temp_boundary_setings.params = {0, 1, flux};
            settings.mesh_settings.boundary_settings_info[2*temp_boundary_setings.axis + temp_boundary_setings.side] = temp_boundary_setings;
            fixedflux_boundaries_element = fixedflux_boundaries_element->NextSiblingElement("FixedFlux");
        }
        TiXmlElement* fixedtemperature_boundaries_element = boundaries_element->FirstChildElement("FixedTemperature");
        while (fixedtemperature_boundaries_element)
        {
            Settings::MeshSettings::BoundarySettings temp_boundary_setings;
            temp_boundary_setings.type = Settings::MeshSettings::BoundarySettings::BoundarySettingTypes::FixedTemperature;
            ParseBoundaryAxis(fixedtemperature_boundaries_element, temp_boundary_setings.axis);
            ParseUnsigned(fixedtemperature_boundaries_element, "side", &temp_boundary_setings.side);
            double temperature;
            ParseScalar(fixedtemperature_boundaries_element, "temperature", &(temperature));
            temp_boundary_setings.params = {1, 0, temperature};
            settings.mesh_settings.boundary_settings_info[2*temp_boundary_setings.axis + temp_boundary_setings.side] = temp_boundary_setings;
            fixedtemperature_boundaries_element = fixedtemperature_boundaries_element->NextSiblingElement("FixedTemperature");
        }
    }
}


void ParseSnapshotInfo(TiXmlElement* snapshot_info_element, Settings& settings)
{
    ParseUnsigned(snapshot_info_element->FirstChildElement("Period"), "frames", &settings.snapshot_settings_info.period_frames);
    TiXmlElement* snapshot_data_element = snapshot_info_element->FirstChildElement("Data");
    ParseString(snapshot_data_element, "fileName", &settings.snapshot_settings_info.snaphot_data_path);
    ParseVector(snapshot_data_element, "stride", &settings.snapshot_settings_info.stride);
    ParseBool(snapshot_data_element, "writeTemperature", &settings.snapshot_settings_info.write_temperature);
    ParseBool(snapshot_data_element, "writeEnthalpy", &settings.snapshot_settings_info.write_enthalpy);
    ParseBool(snapshot_data_element, "writeThermalElasticity", &settings.snapshot_settings_info.write_thermal_elasticity);
    ParseBool(snapshot_data_element, "writeState", &settings.snapshot_settings_info.write_state);
}


void ParseTaskInfo(TiXmlElement* task_info_element, Settings& settings)
{
    ParseScalar(task_info_element, "timeStep", &(settings.task_settings.time_step));
    ParseUnsigned(task_info_element, "numberOfSteps", &(settings.task_settings.number_of_steps));
    TiXmlElement * ini_state_element = task_info_element->FirstChildElement("IniState");

    TiXmlElement * ini_state_per_node_element = ini_state_element->FirstChildElement("PerNode");
    TiXmlElement * ini_state_per_submesh_element = ini_state_element->FirstChildElement("PerSubmesh");
    if (! (!ini_state_per_node_element ^ !ini_state_per_submesh_element))
    {
        std::cerr << "Error: Wrong initial state settings \n";
        std::exit(EXIT_FAILURE);
    }
    else if (ini_state_per_node_element)
    {
        settings.task_settings.initial_state_settings.type = 
            Settings::TaskSettings::InitialStateSettings::InitialStateSettingsType::PerNode;
        ParseString(ini_state_per_node_element, "fileName", &settings.task_settings.initial_state_settings.initial_state_data_file);
    }
    else if (ini_state_per_submesh_element)
    {
        settings.task_settings.initial_state_settings.type = Settings::TaskSettings::InitialStateSettings::InitialStateSettingsType::PerSubmesh;
        TiXmlElement * submesh_element = ini_state_per_submesh_element->FirstChildElement("Submesh");
        while (submesh_element)
        {
            int index;
            Settings::MeshSettings::MediumParams temp_medium_params;
            double temp_temperature;
            ParseScalar(submesh_element, "temperature", &(temp_temperature));
            settings.task_settings.initial_state_settings.initial_temperatures_by_submesh.push_back(temp_temperature);
            submesh_element = submesh_element->NextSiblingElement("Submesh");
        }
    }
}


