####DAG modelling
#dag for main objective of whether occupations with higher 
#blast exposure have worse cognition, mental health, and concussion symptoms
#the primary insight however, is whether cognition is mediated by mental health

#Libraries required
library(dagitty)
library(ggdag)
library(ggplot2)

#Simpler model from which manuscript statistical models are derived
#DAG 1 including RPQ scores
dag1 <- dagify(NCog ~ Occ, 
               MH ~ Occ, NCog ~ MH, MH ~ Age, NCog ~ Age,
               exposure = "Occ",
               outcome = "NCog")

dagitty::is.dagitty(dag1)          
impliedConditionalIndependencies(dag1)  
adjustmentSets(dag1, exposure = "Occ", outcome = "NCog", effect = "total")
adjustmentSets(dag1, exposure = "Occ", outcome = "NCog", effect = "direct")

#checking for mh
adjustmentSets(dag1, exposure = "Occ", outcome = "MH", effect = "total")
adjustmentSets(dag1, exposure = "Occ", outcome = "MH", effect = "direct")

#find mediation paths for the cog model
med_path <- paths(dag1, from = "Occ", to = "NCog")
print(med_path)

#set seed so DAG is the same each time
set.seed(2)    

#draw the DAG with ggdag with smaller text inside nodes
d1 <- ggdag(dag1) + theme_dag() 
d1 <- ggplot(dag1, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_node() + 
  geom_dag_edges() +
  geom_dag_text(size = 3) +
  theme_dag()
d1

#save the DAG
ggsave("fig1.jpg", d1, dpi = 600)


