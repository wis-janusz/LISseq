curl -fsSL https://pyenv.run | bash

export PYENV_ROOT="$HOME/.pyenv"
[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init - bash)"

pyenv install 3.12
pyenv virtualenv 3.12 lisseq
pyenv activate lisseq
pip install -r requirements.txt --disable-pip-version-check
chmod 754 ./run.sh
