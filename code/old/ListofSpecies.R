mutate(curve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                            species == 'Aedes aegypti' ~ '2',
                            species == 'Anthonomus grandis' ~ '3',
                            species == 'Paracoccus marginatus' ~ '4',
                            species == 'Acyrthosiphon pisum' ~ '5',
                            species == 'Aphis gossypii' ~ '6',
                            species == 'Harmonia axyridis' ~ '7',
                            species == 'Tribolium castaneum' ~ '8',
                            species == 'Aedes krombeini' ~ '9',
                            species == 'Bemisia tabaci' ~ '10',
                            species == 'Tetraneura nigriabdominalis' ~ '11',
                            species == 'Stethorus punctillum' ~ '12',
                            species == 'Tetranychus mcdanieli' ~ '13',
                            species == 'Tetranychus urticae' ~ '14',
                            species == 'Clavigralla tomentosicollis' ~ '15',
                            species == 'Planococcus citri' ~ '16',
                            species == 'Muscidifurax zaraptor' ~ '17',
                            species == 'Aphis nasturtii' ~ '18',
                            species == 'Rhopalosiphum maidis' ~ '19',
                            species == 'Anopheles gambiae' ~ '20')) %>%
  
 # need fitting individually:
  
  
  #,
  #  species ==   'Clavigralla tomentosicolli' ~ '8')) %>%
  #  species == 'Paracoccus marginatus' ~ '13',
  