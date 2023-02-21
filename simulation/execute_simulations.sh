# HNUTS no split
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts666 -m e -M achin23@jhu.edu -v seed=666,output=hbps_hnuts666,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts667 -m e -M achin23@jhu.edu -v seed=667,output=hbps_hnuts667,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts668 -m e -M achin23@jhu.edu -v seed=668,output=hbps_hnuts668,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts669 -m e -M achin23@jhu.edu -v seed=669,output=hbps_hnuts669,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts670 -m e -M achin23@jhu.edu -v seed=670,output=hbps_hnuts670,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh


# HNUTS split
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts666 -m e -M achin23@jhu.edu -v seed=666,output=hnuts666,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts667 -m e -M achin23@jhu.edu -v seed=667,output=hnuts667,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts668 -m e -M achin23@jhu.edu -v seed=668,output=hnuts668,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts669 -m e -M achin23@jhu.edu -v seed=669,output=hnuts669,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts670 -m e -M achin23@jhu.edu -v seed=670,output=hnuts670,xml=24t_HNUTS.xml phylogenetic.sh


# BPS
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps666 -m e -M achin23@jhu.edu -v seed=666,output=bps666,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps667 -m e -M achin23@jhu.edu -v seed=667,output=bps667,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps668 -m e -M achin23@jhu.edu -v seed=668,output=bps668,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps669 -m e -M achin23@jhu.edu -v seed=669,output=bps669,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps670 -m e -M achin23@jhu.edu -v seed=670,output=bps670,xml=24t_BPS.xml phylogenetic.sh
