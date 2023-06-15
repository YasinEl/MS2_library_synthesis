using EzXML
using DataFrames
using CSV


directory = "set1" #Directory to all files
file_list = sort(readdir(directory)) #List all files

df = DataFrame(id = String[], name = String[], formula = String[], charge = String[], compartment = String[], reactant_reaction_id = String[], reactant_reaction_name = String[], reactant_reaction_num = Integer[], product_reaction_id = String[], product_reaction_name = String[], product_reaction_num = Integer[], genera = String[], species = String[])

id_and_index = Dict()

for file in file_list
    
	if collect(split(file, "."))[2] != "xml"
		continue
    	else
        
        	path_to_file = directory * "/" * file
	
		try
		
        		readxml(path_to_file)

		catch e
			continue
		end

		doc = readxml(path_to_file)

        	genera = split(file, "_")[1]
        	specie = split(file, ".xml")[1]
        	current_index = 1
        
        	n = root(doc)
        	model_node = firstnode(n)
        	annotation_node = nextnode(model_node)
        	for sub in eachelement(annotation_node)
            	if sub.name == "listOfSpecies"
        
                	for m in eachelement(sub)
                    	#println(m["metaid"])
                    	split_string = split(m["metaid"], "__")
                    	collection = collect(split_string)
                    	id = collect(split(collection[1], "_"))[2]

                        if haskey(id_and_index, id)
                            row = id_and_index[id]

                            genera_so_far = df[row, :genera]
                            species_so_far = df[row, :species]
            
                            last_specie = collect(split(species_so_far, ";"))[end]

                            if last_specie == specie
                                continue
                            else
                                updated_species = species_so_far * ";" * specie
                                updated_genera = genera_so_far * ";" * genera

                                df[row, :genera] = updated_genera
                                df[row, :species] = updated_species
                            end                        
                        else
                        	compartment = collection[3]
                        	name = m["name"]
			
                			try #If there is an error retrieving any of this, skip this metabolite 

                                		m["fbc:charge"]
                    				m["fbc:chemicalFormula"]

                			catch e
                				continue
                			end

                			charge = m["fbc:charge"]
                        	formula = m["fbc:chemicalFormula"]

                        	push!(df, [id, name, formula, charge, compartment, "", "", 0, "", "", 0, genera, specie])
                        	id_and_index[id] = current_index
                        	current_index+=1
                        end
                	end
            	end
    
            	if sub.name == "listOfReactions"
        
                	for r in eachelement(sub) #Go over each reaction
                    	metaid = collect(split(r["metaid"], "_"))[2] #Get reaction ID
                    	r_name = r["name"]
            
                    	for p in eachelement(r) #Go over each element inside until finding reactants and products
                
                       		if p.name == "listOfReactants" #If this is the list of reactants, iterate over each specie that participates
                    
                            		for s in eachelement(p)
                                		compound = collect(split(s["species"], "_"))[2]
                                		row = id_and_index[compound]
                                		reaction_id_values = df[row, :reactant_reaction_id]
                                		reaction_name_values = df[row, :reactant_reaction_name]
                                		reaction_num_values = df[row, :reactant_reaction_num]
                        
                                		updated_reaction_id_values = reaction_id_values * ";" * metaid
                                		updated_reaction_name_values = reaction_name_values * ";" * r_name
                                		updated_reaction_num_values = reaction_num_values + 1
                        
                                		#Update the values
                                		df[row, :reactant_reaction_id] = updated_reaction_id_values
                                		df[row, :reactant_reaction_name] = updated_reaction_name_values
                                		df[row, :reactant_reaction_num] = updated_reaction_num_values
                               
                            		end
                        	end
                
                        	if p.name == "listOfProducts"

                            		for s in eachelement(p)
                                		compound = collect(split(s["species"], "_"))[2]
                                		row = id_and_index[compound]
                                		reaction_id_values = df[row, :product_reaction_id]
                                		reaction_name_values = df[row, :product_reaction_name]
                                		reaction_num_values = df[row, :product_reaction_num]
                        
                                		updated_reaction_id_values = reaction_id_values * ";" * metaid
                                		updated_reaction_name_values = reaction_name_values * ";" * r_name
                                		updated_reaction_num_values = reaction_num_values + 1
                        
                                		#Update the values
                                		df[row, :product_reaction_id] = updated_reaction_id_values
                                		df[row, :product_reaction_name] = updated_reaction_name_values
                                		df[row, :product_reaction_num] = updated_reaction_num_values
                    
                            		end
                        	end
                
			end
            
                end
    
            end
  
        end     
        
    end

message = "Done with file" * " "* file
print(message)
print("\n")
end

CSV.write("metabolite_reconstruction_list.csv", df)
