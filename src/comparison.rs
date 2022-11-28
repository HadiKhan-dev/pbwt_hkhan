pub fn cov(xi: &Vec<f64>,yi: &Vec<f64>) -> f64 {
    let l = xi.len();
    let mut xsum: f64 = 0.0;
    let mut ysum: f64 = 0.0;

    for i in 0..l {
        if xi[i].is_nan() {
            println!("X issue");
        }

        if yi[i].is_nan() {
            println!("Y issue");
        }
        xsum += xi[i];
        ysum += yi[i];
    }

    let xmu = xsum/(l as f64);
    let ymu = ysum/(l as f64);

    let mut covar = 0.0;

    for i in 0..l {
        covar += (xi[i]-xmu)*(yi[i]-ymu);
    }
    
    return covar/((l as f64)-1.0);
}

pub fn corr(xi: &Vec<f64>,yi: &Vec<f64>) -> f64 {
    return cov(&xi,&yi)/(cov(&xi,&xi)*cov(&yi,&yi)).sqrt();
}


pub fn compare_results(true_values: &Vec<Vec<u8>>, imputed_values: &Vec<Vec<u8>>, allele_freqs: &Vec<f64>, buckets : &Vec<f64>) -> () {
    let mut bucketed: Vec<Vec<[f64;2]>> = vec![Vec::new(); buckets.len()];
    let N = true_values[0].len();
    let M = true_values.len();


    for j in 0..N {
        let freq = allele_freqs[j];
        let position = buckets.binary_search_by(|v| {
            v.partial_cmp(&freq).expect("Couldn't compare values")
        });
        let mut loc;
        match position {
            Ok(i) => {
                loc = i;
            }
            Err(i) => {
                loc = i;
            }
        }

        for i in 0..(M/2) {
            let true_sum = true_values[2*i][j]+true_values[2*i+1][j];
            let imputed_sum = imputed_values[2*i][j]+imputed_values[2*i+1][j];

            bucketed[loc].push([true_sum as f64, imputed_sum as f64]);
        }
        
    }

    for j in 0..bucketed.len() {
        let name = format!("{}",buckets[j]);

        let size = bucketed[j].len();
        let mut trues:  Vec<f64> = Vec::with_capacity(size);
        let mut imputes: Vec<f64> = Vec::with_capacity(size);

        let mut ct_zero_true: f64 = 0.0;
        let mut ct_zero_impute: f64 = 0.0;

        for i in 0..size {
            let true_add = bucketed[j][i][0];
            let imp_add = bucketed[j][i][1];

            if true_add == 0.0 {
                ct_zero_true += 1.0;
            }
            if imp_add == 0.0 {
                ct_zero_impute += 1.0;
            }

            trues.push(bucketed[j][i][0]);
            imputes.push(bucketed[j][i][1]);
        }

        let corr = corr(&trues,&imputes);

        println!("Bucket: {}, Corr: {:.4?}, Len: {}, Bias: {} {} ", name,corr,trues.len()/5,(ct_zero_true as f64)/(trues.len() as f64),(ct_zero_impute as f64)/(trues.len() as f64));

    }


}