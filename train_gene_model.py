# train a miniature character-level synteny model

out_dir = '/media/mirrored-hdd/shorsfield/jobs/pangenome_LLM/models/gene_char_geneID_2k_iter_A'
eval_interval = 100 # keep frequent because we'll overfit
eval_iters = 100
log_interval = 10 # don't print too too often

# we expect to overfit on this small dataset, so only save when val improves
always_save_checkpoint = False

wandb_log = False # override via command line if you like
wandb_project = 'gene_char'
wandb_run_name = 'mini-gpt'

data_dir = '/media/mirrored-hdd/shorsfield/jobs/pangenome_LLM/training_data/gene_char_geneID'
dataset = 'gene_char'
gradient_accumulation_steps = 1
batch_size = 16
block_size = 1024 # context of up to N previous characters, 81K characters at a time with current settings

# based on panGPT and gpt2
n_layer = 12
n_head = 12
n_embd = 768
dropout = 0.2

learning_rate = 1e-4
max_iters = 2000 # equal to N epochs of 85 million token dataset
lr_decay_iters = 2000 # make equal to max_iters usually
min_lr = 1e-5 # learning_rate / 10 usually
#beta2 = 0.99 # make a bit bigger because number of tokens per iter is small

warmup_iters = 100 # not super necessary potentially

# weight decay
weight_decay = 1e-5

# on macbook also add
# device = 'cpu'  # run on cpu only
# compile = False # do not torch compile the model
