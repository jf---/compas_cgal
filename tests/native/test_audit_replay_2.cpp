void audit_replay_transaction_state_gate();
void audit_replay_trial_depletion_gate();
void audit_replay_identity_lineage_gate();
void audit_replay_finalization_gate();

void audit_replay_gate()
{
    audit_replay_transaction_state_gate();
    audit_replay_trial_depletion_gate();
    audit_replay_identity_lineage_gate();
    audit_replay_finalization_gate();
}
