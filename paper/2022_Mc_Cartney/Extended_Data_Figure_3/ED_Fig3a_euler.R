pdf(file="0904.self.svs.e.pdf", useDingbats = FALSE, width = 4, height = 4)

library(eulerr)

s2 <- c("CLR" = 2037,
        "HIFI" = 66,
        "ONT" = 3586,
        
        "CLR&HIFI" = 5,
        "CLR&ONT" = 12,
        "HIFI&ONT" = 10,
        
        "CLR&HIFI&ONT" = 17
)

#plot(venn(s2))
plot(euler(s2), quantities = TRUE, shape = "ellipse")

dev.off()
