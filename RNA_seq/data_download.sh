# create a new folder anywhere, e.g. Desktop (use you own login name to substitute your_user_name)
mkdir -p '/mnt/c/Users/your_user_name/Desktop/unix'

# link the folder to linux system
ln -s  '/mnt/c/Users/your_user_name/Desktop/unix' ~/unix

# check
cd
ls -lh

# success if you see this line, from now on, you could manipulate files in both systems
# lrwxrwxrwx 1 USERNAME USERNAME  DATE TIME unix -> /mnt/c/Users/your_user_name/Desktop/unix

# Please download and install firehose_get first!
# Then,
firehose_get -tasks rnaseq_pre stddata latest blca

# enter the inner directory
cd stddata__2016_07_15/BLCA/20160715/

# check files
ls -lh

# delete md5, aux and mage-tab files
rm -f *md5 *aux* *mage-tab*
ls -lh

# decompression
tar zxvf gdac.broadinstitute.org_BLCA.mRNAseq_Preprocess.Level_3.2016071500.0.0.tar.gz
cd gdac.broadinstitute.org_BLCA.mRNAseq_Preprocess.Level_3.2016071500.0.0/
ls -lh

# we only need BLCA.uncv2.mRNAseq_raw_counts.txt, see details below
mkdir rnaseq_test
cp BLCA.uncv2.mRNAseq_raw_counts.txt ./rnaseq_test/
cd rnaseq_test/
ls -lh
