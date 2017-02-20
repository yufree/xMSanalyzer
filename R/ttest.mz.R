ttest.mz <-
function(x, y)
				                                                        {
				                                                                
				                                                                
				                                                               
				                                                                        ttest_res=t.test(t(x), as.matrix(y),paired=T)
													ttest_res=abs(ttest_res$p.value)
				                                                                
				                                                                return(ttest_res)
				                                                        }
