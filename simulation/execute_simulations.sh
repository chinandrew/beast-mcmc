# HNUTS no split
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts666 -m e -M achin23@jhu.edu -v seed=666,output=hbps_hnuts666,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts667 -m e -M achin23@jhu.edu -v seed=667,output=hbps_hnuts667,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts668 -m e -M achin23@jhu.edu -v seed=668,output=hbps_hnuts668,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts669 -m e -M achin23@jhu.edu -v seed=669,output=hbps_hnuts669,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts670 -m e -M achin23@jhu.edu -v seed=670,output=hbps_hnuts670,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts671 -m e -M achin23@jhu.edu -v seed=671,output=hbps_hnuts671,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts672 -m e -M achin23@jhu.edu -v seed=672,output=hbps_hnuts672,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts673 -m e -M achin23@jhu.edu -v seed=673,output=hbps_hnuts673,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts674 -m e -M achin23@jhu.edu -v seed=674,output=hbps_hnuts674,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hbps_hnuts675 -m e -M achin23@jhu.edu -v seed=675,output=hbps_hnuts675,xml=24t_HBPS_NUTS_0-001.xml phylogenetic.sh


# HNUTS split
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts666 -m e -M achin23@jhu.edu -v seed=666,output=hnuts666,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts667 -m e -M achin23@jhu.edu -v seed=667,output=hnuts667,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts668 -m e -M achin23@jhu.edu -v seed=668,output=hnuts668,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts669 -m e -M achin23@jhu.edu -v seed=669,output=hnuts669,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts670 -m e -M achin23@jhu.edu -v seed=670,output=hnuts670,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts671 -m e -M achin23@jhu.edu -v seed=671,output=hnuts671,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts672 -m e -M achin23@jhu.edu -v seed=672,output=hnuts672,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts673 -m e -M achin23@jhu.edu -v seed=673,output=hnuts673,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts674 -m e -M achin23@jhu.edu -v seed=674,output=hnuts674,xml=24t_HNUTS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_hnuts675 -m e -M achin23@jhu.edu -v seed=675,output=hnuts675,xml=24t_HNUTS.xml phylogenetic.sh


# BPS
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps666 -m e -M achin23@jhu.edu -v seed=666,output=bps666,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps667 -m e -M achin23@jhu.edu -v seed=667,output=bps667,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps668 -m e -M achin23@jhu.edu -v seed=668,output=bps668,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps669 -m e -M achin23@jhu.edu -v seed=669,output=bps669,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps670 -m e -M achin23@jhu.edu -v seed=670,output=bps670,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps671 -m e -M achin23@jhu.edu -v seed=671,output=bps671,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps672 -m e -M achin23@jhu.edu -v seed=672,output=bps672,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps673 -m e -M achin23@jhu.edu -v seed=673,output=bps673,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps674 -m e -M achin23@jhu.edu -v seed=674,output=bps674,xml=24t_BPS.xml phylogenetic.sh
qsub -q "shared.q@compute-124"  -l mem_free="16G",h_vmem="16G" -N p_bps675 -m e -M achin23@jhu.edu -v seed=675,output=bps675,xml=24t_BPS.xml phylogenetic.sh
