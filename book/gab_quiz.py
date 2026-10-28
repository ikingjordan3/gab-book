"""Self-check quizzes in the course colours. Used by every module page:

    import sys; sys.path.insert(0, "..")
    from gab_quiz import quiz
    quiz("quizzes/<module>.json")
"""
from jupyterquiz import display_quiz

NAVY, GREEN, RED = "#3B5B92", "#2E7D5B", "#B03A2E"
COLORS = {
    "--jq-multiple-choice-bg": NAVY,
    "--jq-many-choice-bg": NAVY,
    "--jq-numeric-bg": NAVY,
    "--jq-string-bg": NAVY,
    "--jq-correct-color": GREEN,
    "--jq-incorrect-color": RED,
}


def quiz(path):
    """Show a re-runnable quiz: answers are shuffled, feedback explains every option, nothing is recorded."""
    display_quiz(path, shuffle_answers=True, colors=COLORS)
