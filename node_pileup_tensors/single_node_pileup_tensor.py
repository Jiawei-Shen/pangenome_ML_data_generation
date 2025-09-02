def display_pileup_c1_only(variant_key: str,
                           v_pos: int,
                           v_type: str,
                           ch1: np.ndarray,
                           reads_added: int):
    """
    Print only Channel 1 (bases) for ref and the first `reads_added` reads.
    ch1 shape: (1 + TENSOR_MAX_READ_ROWS, TENSOR_WINDOW_SIZE)
    """
    # compute variant column in the window
    win_center = v_pos + 1 if v_type == 'I' else v_pos
    win_start  = calculate_window_start(win_center, TENSOR_WINDOW_SIZE)
    v_idx      = v_pos - win_start

    def _row_to_bases_str(int_row):
        return ''.join(INDEX_TO_BASE_FOR_VIEW.get(int(x), '?') for x in int_row)

    print(f"\n--- Variant: {variant_key} ---")

    # Ref row
    ref_c1 = _row_to_bases_str(ch1[0, :])
    print(f"  Ref C1: {ref_c1}")
    # caret marker under the variant column
    marker = ''.join('^' if i == v_idx else ' ' for i in range(TENSOR_WINDOW_SIZE))
    print(f"       ^: {marker}")

    # Read rows (only C1)
    for r in range(1, 1 + reads_added):
        row_c1 = _row_to_bases_str(ch1[r, :])
        print(f"  R{r:03d} C1: {row_c1}")
